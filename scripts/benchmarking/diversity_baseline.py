"""Designs made the current way, measured per strain and per host.

Phase 3 of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.
The question is whether a primer set designed the way the tool designs one
today loses coverage on strains it was not designed against, or specificity
against hosts it was not designed against. This script builds the designs and
measures them. It decides nothing, and it changes no product code.

WHAT IS BUILT. Every design is built from `count-kmers` onward in its own data
directory, because the prepared example under `examples/wolbachia_pool_design`
cannot be re-optimized (its position indexes predate reference-digest
recording, and `optimize` refuses instead of guessing).

  D1  reference strain only, one host            today's usual run
  D2  every strain pooled in fg_genomes, one host today's multi-genome run
  C1  D2 with one strain left out                 evaluated on the strain left out
  D3  reference strain only, every allowed host   pooled hosts
  C2  random size-matched panels drawn from the candidate pool of a design,
      one per seed; no design is built for these

Steps 1 to 3 do not read the panel size, so one candidate pool serves every
requested size and only `optimize` is repeated. Its outputs are moved to
`sets/nNN/` inside the data directory after each run.

WHAT IS MEASURED. Set 0 of every delivered design, and every C2 panel, goes
through `reference_panel_evaluation.evaluate_reference_panel`, the Phase 2
function, against every strain and every allowed host. A reference that took
part in the design is read through the design's own prefix (its position index
holds the whole candidate pool); any other reference is read through a prefix
under `eval/`, where only a k-mer table exists and positions come from scanning
the FASTA. The record says which route each figure used.

UNKNOWN IS NOT ZERO. A figure that could not be measured is carried as the
`Measurement` the Phase 2 function returned, with its reason. A design that
fails is a result: its `design_failure.json` is recorded and the run moves on.
A summary over seeds with one unavailable member is unavailable.

`num_primers` IS A REQUEST. The delivered panel may be smaller. Requested and
delivered sizes are recorded separately, and C2 panels match the DELIVERED
size.

C2 PANELS ARE NOT DIMER-SCREENED. They are uniform draws from `step3_df.csv`.
They answer "what does any panel from this pool do", not "what does any
orderable panel do".

MEMORY. Each pipeline step and each evaluation batch runs as a child process
under a watchdog. The watchdog samples the resident set size of the child's
process tree (from `ps`, summed over the tree, so pages shared between workers
are counted more than once and the figure is an upper bound; the largest single
process is recorded beside it) and the
system-wide free percentage (from `memory_pressure`). It stops the child when
either limit is passed. A stopped step is recorded with its peak and is not
retried unless `--retry` is given. No step starts while the free percentage is
below the start threshold. Steps run strictly one at a time.

THE LARGE HOST. A host longer than `--max-host-bp` is refused as a design
background, so the 3.3 Gb human reference cannot enter a design by accident.

A host longer than `COUNTS_ONLY_ABOVE_BP` is handled differently again: it is
admitted as an EVALUATION REFERENCE on the counts route and as nothing else.
The rule is the length, not a flag, so no combination of flags reaches the
other routes:

  - `load_panel` puts it in `counts_only_hosts`, never in `hosts`, whichever
    flag named it. `plan_designs` reads `hosts` alone, so no design can take it
    as a background and no `filter` run ever sees it.
  - `design_params` refuses to write parameters that name one, which is the
    second check on the same fact.
  - `reference_plan` sets `scan=False` for one unconditionally, ignoring
    `--c2-scan-offdesign-hosts`, so its positions are never looked for.

Its sites and densities are then measured exactly and its coverage and gap
figures are unavailable with the reason the Phase 2 function gives.

THE COUNTS-ONLY PASS. `--counts-only-pass` adds those references to an existing
results file: it reads `results.json`, measures every delivered set 0 and every
C2 panel already recorded there against the counts-only hosts, and writes the
figures back beside the ones already present. It builds no design and rebuilds
no pool. A design or a step record it did not produce is carried through
untouched, and it refuses a results file that holds no design.

RESUMING. A step whose record says it finished, under the same parameters, and
whose artifacts are present is reused. A data directory built under different
parameters is refused instead of overwritten.

Usage (a k-mer counter must be on PATH):

    python scripts/benchmarking/diversity_baseline.py \
        --manifest tests/validation/genomes/diversity_panel.json \
        --workdir tests/validation/genomes/diversity_baseline

    # one group of designs, or fewer sizes
    python scripts/benchmarking/diversity_baseline.py --only D1 --sizes 6
"""

import argparse
import csv
import hashlib
import json
import os
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

DEFAULT_MANIFEST = os.path.join("tests", "validation", "genomes", "diversity_panel.json")
DEFAULT_WORKDIR = os.path.join("tests", "validation", "genomes", "diversity_baseline")
DEFAULT_HOSTS = ("lactobacillus", "drosophila")
DEFAULT_SIZES = (6, 12, 24)
GROUPS = ("D1", "D2", "C1", "C2", "D3")

K = 12

# A reference longer than this is evaluated by k-mer counts and by nothing
# else, whatever flag named it. Reaching positions in a reference means holding
# it in memory (`query_scan` is measured at 2.3 GB on hg38) and a `filter` run
# against hg38 peaks at about 8.5 GB, neither of which belongs in this stage.
# It is a constant rather than an option so that no flag combination can send
# such a reference down the position route or into a design background.
COUNTS_ONLY_ABOVE_BP = 1_000_000_000

# The settings of examples/wolbachia_pool_design/params.template.json, which is
# the pool the project's recorded figures were measured on. Only the genomes,
# the prefixes, the data directory and the panel size vary between designs.
DESIGN_PARAMS = {
    "schema_version": 2,
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": K,
    "max_k": K,
    "gc_tolerance": 0.15,
    "gc_min": 0.2,
    "gc_max": 0.5,
    "dmso_percent": 0.0,
    "betaine_m": 0.0,
    "trehalose_m": 0.0,
    "na_conc": 50.0,
    "mg_conc": 10.0,
    "min_fg_freq": 1e-06,
    "max_bg_freq": 5e-06,
    "max_gini": 0.7,
    "max_primer": 2000,
    "min_tm": 25.0,
    "max_tm": 55.0,
    "max_dimer_bp": 3,
    "max_self_dimer_bp": 4,
    "optimization_method": "hybrid",
    "iterations": 8,
    "max_sets": 5,
    "cpus": 4,
    "verbose": True,
    "fg_circular": True,
    "bg_circular": False,
}

# The template states the reference's GC content. It is kept for a design
# whose foreground is the reference alone and left to the pipeline otherwise.
REFERENCE_GENOME_GC = 0.3523

# The project's recorded figure for the Wolbachia pool, and where it came from.
RECORDED_REFERENCE = {
    "selectivity_density": 60.112,
    "coverage": 0.6535,
    "panel_size": 12,
    "path": "plan-pool with a selectivity-density floor of 60, occupancy-weighted coverage",
    "pool": "examples/wolbachia_pool_design/work/step3_df.csv (wMel against Drosophila)",
    "source": "docs/validation/frontier_refill_2026-09-17.md, "
    "docs/validation/no_search_headroom_on_this_pool_2026-09-18.md",
}

OPTIMIZE_SEED = 0

STATUS_OK = "ok"
STATUS_FAILED = "failed"
STATUS_STOPPED = "stopped"
STATUS_REUSED = "reused"

EXIT_HALTED = 3


class RunHalted(Exception):
    """The run must stop: memory stayed low, or a step passed the halt peak."""


# ----------------------------------------------------------------------
# Memory: what the system has, what a child tree holds, and the watchdog
# ----------------------------------------------------------------------


def free_percent():
    """System-wide free memory percentage from `memory_pressure`, or None.

    None means the figure could not be read. It is never treated as plenty:
    the callers refuse to start, and the watchdog reports it as unread.
    """
    try:
        out = subprocess.run(
            ["memory_pressure"], capture_output=True, text=True, timeout=20, check=False
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    for line in reversed(out.splitlines()):
        if "free percentage" in line:
            try:
                return int(line.rsplit(":", 1)[1].strip().rstrip("%"))
            except ValueError:
                return None
    return None


def tree_rss_bytes(root_pid):
    """Resident set size of a process tree, or None when it could not be read.

    Returns `(total, largest, processes)`: the sum over the root and its
    descendants, the largest single process among them, and how many there
    were. Read from `ps`. Shared pages are counted once per process, so for a
    step that forks workers the total is an upper bound on what the tree
    holds; the largest single process says how much of it is one process.
    """
    try:
        out = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    children = {}
    rss = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            pid, ppid, kb = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            continue
        children.setdefault(ppid, []).append(pid)
        rss[pid] = kb
    if root_pid not in rss:
        return None
    total = 0
    largest = 0
    stack = [root_pid]
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += rss.get(pid, 0)
        largest = max(largest, rss.get(pid, 0))
        stack.extend(children.get(pid, []))
    return total * 1024, largest * 1024, len(seen)


def wait_for_memory(limits, label):
    """Block until enough memory is free to start a step, or halt the run.

    Free disk is checked first and is not waited for: a step that fills the
    disk does not recover by waiting.
    """
    free_disk = shutil.disk_usage(limits["disk_path"]).free
    if free_disk < limits["min_free_disk_bytes"]:
        raise RunHalted(
            f"{free_disk / 2**30:.1f} GiB of disk is free before {label}, below "
            f"{limits['min_free_disk_bytes'] / 2**30:.1f} GiB"
        )
    deadline = time.time() + limits["start_wait_seconds"]
    while True:
        free = free_percent()
        if free is not None and free >= limits["start_free_percent"]:
            return free
        if time.time() >= deadline:
            raise RunHalted(
                f"free memory was {free} percent before {label} and stayed below "
                f"{limits['start_free_percent']} for {limits['start_wait_seconds']} s"
            )
        print(f"  waiting: {free} percent free before {label}", flush=True)
        time.sleep(30)


def _downsample(samples, keep):
    """At most `keep` samples, evenly spaced, with the last one always kept.

    The peak figures are computed over every sample; this is only what gets
    written down, so a long step does not put thousands of rows in the record.
    """
    if len(samples) <= keep:
        return samples
    step = len(samples) / float(keep)
    picked = [samples[int(i * step)] for i in range(keep)]
    if picked[-1] is not samples[-1]:
        picked[-1] = samples[-1]
    return picked


def run_watched(cmd, log_path, limits, label):
    """Run one child under the watchdog and return what it cost.

    The child gets its own session so the whole tree can be signalled. The
    record carries the peak tree RSS, the wall time and the free percentage
    before, lowest during and after. A child stopped by the watchdog has
    status `stopped` and the reason; it is not run again here.

    The peak is also split by how many processes the tree held at each sample.
    A step that invokes an external counter spends its first phase as the
    `neoswga` process alone and its second with the counter as a child, and the
    two phases can peak for completely different reasons: on hg38 the first
    phase peaked at 6.18 GiB holding two copies of the reference, before the
    counter had started at all. One figure over the whole step hides which
    phase the memory belongs to, so `peak_alone_bytes` and `peak_with_children_bytes`
    report them apart. A downsampled trace is kept beside them so the shape can
    be read rather than inferred from three numbers.
    """
    free_before = wait_for_memory(limits, label)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    started = time.time()
    peak = None
    at_peak = None
    peak_single = None
    most_processes = None
    # Peaks by phase. None means that phase was never sampled, which is not a
    # peak of zero: a step whose child starts before the first sample has no
    # observed alone phase.
    peak_alone = None
    peak_with_children = None
    trace = []
    lowest_free = free_before
    stop_reason = None
    with open(log_path, "w") as log:
        child = subprocess.Popen(
            cmd, cwd=REPO_ROOT, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
        )
        while True:
            try:
                child.wait(timeout=limits["sample_seconds"])
                break
            except subprocess.TimeoutExpired:
                pass
            sample = tree_rss_bytes(child.pid)
            rss = None
            if sample is not None:
                rss, largest, processes = sample
                if peak is None or rss > peak:
                    peak = rss
                    at_peak = {"largest_process_rss_bytes": largest, "processes": processes}
                peak_single = largest if peak_single is None else max(peak_single, largest)
                most_processes = (
                    processes if most_processes is None else max(most_processes, processes)
                )
                if processes <= 1:
                    peak_alone = rss if peak_alone is None else max(peak_alone, rss)
                else:
                    peak_with_children = (
                        rss if peak_with_children is None else max(peak_with_children, rss)
                    )
                trace.append(
                    {
                        "at_seconds": round(time.time() - started, 1),
                        "tree_rss_bytes": rss,
                        "processes": processes,
                    }
                )
            free = free_percent()
            if free is not None:
                lowest_free = free if lowest_free is None else min(lowest_free, free)
            if rss is not None and rss > limits["rss_limit_bytes"]:
                stop_reason = (
                    f"tree RSS {rss / 2**30:.2f} GiB passed the limit of "
                    f"{limits['rss_limit_bytes'] / 2**30:.2f} GiB"
                )
            elif free is not None and free < limits["stop_free_percent"]:
                stop_reason = (
                    f"system free memory fell to {free} percent, below "
                    f"{limits['stop_free_percent']}"
                )
            if stop_reason:
                _stop_tree(child)
                break
    wall = time.time() - started
    if stop_reason:
        status = STATUS_STOPPED
    elif child.returncode == 0:
        status = STATUS_OK
    else:
        status = STATUS_FAILED
    return {
        "status": status,
        "returncode": child.returncode,
        "stop_reason": stop_reason,
        "wall_seconds": round(wall, 2),
        # None when the child ended before the first sample, which is a
        # peak that was not observed and not a peak of zero.
        "peak_tree_rss_bytes": peak,
        # The tree at the sample where the total peaked, and the largest any
        # one process reached at any sample. Together they say whether a peak
        # is one process or several holding a copy each.
        "tree_at_peak": at_peak,
        "peak_single_process_rss_bytes": peak_single,
        "most_processes": most_processes,
        # Split by phase. `alone` is the neoswga process before any external
        # counter has been forked; `with_children` is the tree once one has.
        # None means that phase was never sampled.
        "peak_alone_bytes": peak_alone,
        "peak_with_children_bytes": peak_with_children,
        "trace": _downsample(trace, 120),
        "sample_seconds": limits["sample_seconds"],
        "free_percent_before": free_before,
        "free_percent_lowest": lowest_free,
        "free_percent_after": free_percent(),
        "log": os.path.relpath(log_path, REPO_ROOT),
        "command": [os.path.relpath(c, REPO_ROOT) if os.path.isabs(c) else c for c in cmd],
    }


def _stop_tree(child):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(child.pid, sig)
        except (ProcessLookupError, PermissionError):
            break
        try:
            child.wait(timeout=10)
            break
        except subprocess.TimeoutExpired:
            continue


def check_halt(record, limits, label):
    """Halt the run when a step was stopped or its peak passed the halt figure."""
    if record["status"] == STATUS_STOPPED:
        raise RunHalted(f"{label} was stopped by the watchdog: {record['stop_reason']}")
    peak = record.get("peak_tree_rss_bytes")
    if peak is not None and peak > limits["halt_peak_bytes"]:
        raise RunHalted(
            f"{label} peaked at {peak / 2**30:.2f} GiB, above the halt figure of "
            f"{limits['halt_peak_bytes'] / 2**30:.2f} GiB"
        )


# ----------------------------------------------------------------------
# The manifest
# ----------------------------------------------------------------------


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_panel(manifest_path, host_keys, max_host_bp, counts_only_keys=()):
    """Targets and allowed hosts from the manifest, each verified by SHA-256.

    A host that is not named is not opened at all. A named host longer than
    `COUNTS_ONLY_ABOVE_BP` becomes a counts-only evaluation reference whichever
    list named it, and is returned apart from the design hosts; the length
    decides this, not the flag, so a host that large cannot reach a design
    background or a position scan. Among the rest, one longer than
    `max_host_bp` is refused.
    """
    with open(manifest_path) as handle:
        manifest = json.load(handle)
    entries = {entry["key"]: entry for entry in manifest["panel"]}
    targets = sorted(
        (e for e in entries.values() if e.get("role") == "target"), key=lambda e: e["key"]
    )
    if not targets:
        raise ValueError(f"{manifest_path} lists no entry with role 'target'.")
    references = [e for e in targets if e.get("reference")]
    if len(references) != 1:
        raise ValueError(
            f"{manifest_path} must mark exactly one target as the reference, "
            f"found {len(references)}."
        )
    hosts = []
    counts_only = []
    named = list(dict.fromkeys(list(host_keys) + list(counts_only_keys)))
    for key in named:
        entry = entries.get(key)
        if entry is None or entry.get("role") != "host":
            raise ValueError(f"{key!r} is not a host in {manifest_path}.")
        if int(entry["length"]) > COUNTS_ONLY_ABOVE_BP:
            # Whichever list named it. A reference this long is an evaluation
            # reference on the counts route and nothing else.
            counts_only.append(entry)
            continue
        if key in counts_only_keys and key not in host_keys:
            counts_only.append(entry)
            continue
        if int(entry["length"]) > max_host_bp:
            raise ValueError(
                f"host {key!r} is {int(entry['length']):,} bp, above --max-host-bp "
                f"{max_host_bp:,}. A host of that size is a separate stage with its "
                f"own memory plan and is not run from here by default."
            )
        hosts.append(entry)
    skipped = sorted(
        e["key"] for e in entries.values() if e.get("role") == "host" and e["key"] not in named
    )

    panel = {}
    for entry in targets + hosts + counts_only:
        path = entry["path"]
        if not os.path.isabs(path):
            path = os.path.join(REPO_ROOT, path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{entry['key']}: {path} is not on disk.")
        measured = sha256_of(path)
        if measured != entry["sha256"]:
            raise ValueError(
                f"{entry['key']}: {path} has SHA-256 {measured}, the manifest records "
                f"{entry['sha256']}. The genome on disk is not the one the panel describes."
            )
        panel[entry["key"]] = {
            "key": entry["key"],
            "path": path,
            "role": entry["role"],
            "length": int(entry["length"]),
            "records": entry.get("records"),
            "reference": bool(entry.get("reference", False)),
            "sha256": entry["sha256"],
            "counts_only": entry in counts_only,
        }
    return {
        "genomes": panel,
        "targets": [e["key"] for e in targets],
        "hosts": [e["key"] for e in hosts],
        # Deliberately not merged into "hosts": that list is what the designs
        # are planned from, and keeping these out of it is the guarantee that
        # no design takes one as a background.
        "counts_only_hosts": [e["key"] for e in counts_only],
        "reference": references[0]["key"],
        "hosts_not_run": skipped,
    }


# ----------------------------------------------------------------------
# The designs
# ----------------------------------------------------------------------


def plan_designs(panel, c1_host):
    """Every design of this stage, smallest first within each group."""
    reference = panel["reference"]
    targets = panel["targets"]
    hosts = sorted(panel["hosts"], key=lambda key: panel["genomes"][key]["length"])
    designs = []
    for host in hosts:
        designs.append(
            {
                "id": f"D1__{reference}__vs__{host}",
                "group": "D1",
                "fg": [reference],
                "bg": [host],
                "held_out": None,
                "c2": True,
            }
        )
    for host in hosts:
        designs.append(
            {
                "id": f"D2__pooled__vs__{host}",
                "group": "D2",
                "fg": list(targets),
                "bg": [host],
                "held_out": None,
                "c2": True,
            }
        )
    for held in targets:
        designs.append(
            {
                "id": f"C1__without_{held}__vs__{c1_host}",
                "group": "C1",
                "fg": [key for key in targets if key != held],
                "bg": [c1_host],
                "held_out": held,
                "c2": False,
            }
        )
    if len(hosts) > 1:
        designs.append(
            {
                "id": f"D3__{reference}__vs__" + "+".join(hosts),
                "group": "D3",
                "fg": [reference],
                "bg": list(hosts),
                "held_out": None,
                "c2": True,
            }
        )
    return designs


def design_params(design, panel, data_dir, cpus=None):
    genomes = panel["genomes"]
    # The second check on the one fact. `plan_designs` cannot produce such a
    # design, so reaching this means a caller built one by hand; parameters
    # that name a counts-only reference would be read by `filter`.
    forbidden = set(panel.get("counts_only_hosts", ())) & (set(design["fg"]) | set(design["bg"]))
    if forbidden:
        raise ValueError(
            f"design {design['id']!r} names {sorted(forbidden)} as a design genome. "
            f"A reference longer than {COUNTS_ONLY_ABOVE_BP:,} bp is an evaluation "
            f"reference on the counts route only; it is never counted, filtered or "
            f"optimized against here."
        )
    params = dict(DESIGN_PARAMS)
    if cpus is not None:
        params["cpus"] = int(cpus)
    params["data_dir"] = data_dir
    for side in ("fg", "bg"):
        keys = design[side]
        params[f"{side}_genomes"] = [genomes[key]["path"] for key in keys]
        params[f"{side}_prefixes"] = [os.path.join(data_dir, key) for key in keys]
        params[f"{side}_seq_lengths"] = [genomes[key]["length"] for key in keys]
    if design["fg"] == [panel["reference"]]:
        params["genome_gc"] = REFERENCE_GENOME_GC
    params["num_primers"] = DEFAULT_SIZES[0]
    params["target_set_size"] = DEFAULT_SIZES[0]
    return params


def params_digest(params):
    """A digest of the parameters that decide the candidate pool.

    The panel size is left out: steps 1 to 3 do not read it. So is `cpus`,
    which is how many workers a step uses and is recorded on each step
    instead; `--cpus-check` measures whether it changes the pool.
    """
    left_out = ("num_primers", "target_set_size", "cpus")
    pool = {k: v for k, v in params.items() if k not in left_out}
    return hashlib.sha256(json.dumps(pool, sort_keys=True).encode()).hexdigest()


def cli_command(step, params_path, *extra):
    return [sys.executable, "-m", "neoswga.cli_unified", step, "-j", params_path, *extra]


def load_state(data_dir):
    path = os.path.join(data_dir, "baseline_state.json")
    if os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    return {"steps": {}}


def save_state(data_dir, state):
    path = os.path.join(data_dir, "baseline_state.json")
    with open(path + ".tmp", "w") as handle:
        json.dump(state, handle, indent=2)
        handle.write("\n")
    os.replace(path + ".tmp", path)


def pool_artifacts_present(step, params, data_dir):
    """Whether the files a finished step leaves behind are all there."""
    from neoswga.core import kmer_tables
    from neoswga.core.position_index import index_path

    prefixes = params["fg_prefixes"] + params.get("bg_prefixes", [])
    if step == "count-kmers":
        return all(kmer_tables.table_exists(prefix, K) for prefix in prefixes)
    if step == "filter":
        return os.path.exists(os.path.join(data_dir, "step2_df.csv")) and all(
            os.path.exists(index_path(prefix, K)) for prefix in prefixes
        )
    if step == "prepare-candidates":
        return os.path.exists(os.path.join(data_dir, "step3_df.csv"))
    raise ValueError(step)


def run_step(name, cmd, data_dir, state, limits, retry, complete, label, cpus=None):
    """Run one step unless a finished record and its artifacts already exist.

    `cpus` is what the step's params file asks for. It is written on the
    record of a step that runs; a reused record keeps the figure it ran under.
    """
    previous = state["steps"].get(name)
    if previous is not None:
        if previous["status"] == STATUS_OK and complete():
            return dict(previous, reused=True)
        if previous["status"] in (STATUS_STOPPED, STATUS_FAILED) and not retry:
            return dict(previous, reused=True)
    print(f"  {label}: {name}", flush=True)
    record = run_watched(cmd, os.path.join(data_dir, "logs", f"{name}.log"), limits, label)
    record["reused"] = False
    record["cpus"] = cpus
    if record["status"] == STATUS_OK and not complete():
        record["status"] = STATUS_FAILED
        record["stop_reason"] = "the step exited 0 but its artifacts are not all present"
    state["steps"][name] = record
    save_state(data_dir, state)
    peak = record["peak_tree_rss_bytes"]
    peak_text = "not sampled" if peak is None else f"{peak / 2**20:.0f} MiB"
    print(
        f"    {record['status']} in {record['wall_seconds']} s, peak {peak_text}, "
        f"free {record['free_percent_before']} -> {record['free_percent_after']} percent",
        flush=True,
    )
    check_halt(record, limits, f"{label} {name}")
    return record


def build_pool(design, panel, workdir, limits, retry, cpus, data_dir=None):
    """Steps 1 to 3 for one design, in its own data directory."""
    data_dir = data_dir or os.path.join(workdir, "pools", design["id"])
    os.makedirs(data_dir, exist_ok=True)
    params = design_params(design, panel, data_dir, cpus)
    digest = params_digest(params)
    state = load_state(data_dir)
    if state.get("params_digest") not in (None, digest):
        raise ValueError(
            f"{data_dir} was built under different parameters. Remove it, or use "
            f"another --workdir; it is not overwritten."
        )
    state["params_digest"] = digest
    params_path = os.path.join(data_dir, "params.json")
    with open(params_path, "w") as handle:
        json.dump(params, handle, indent=2)
        handle.write("\n")
    save_state(data_dir, state)

    steps = {}
    usable = True
    for step in ("count-kmers", "filter", "prepare-candidates"):
        if not usable:
            steps[step] = {"status": "not_run", "reason": "an earlier step did not finish"}
            continue
        record = run_step(
            step,
            cli_command(step, params_path),
            data_dir,
            state,
            limits,
            retry,
            lambda step=step: pool_artifacts_present(step, params, data_dir),
            design["id"],
            cpus=params["cpus"],
        )
        steps[step] = record
        usable = record["status"] == STATUS_OK
    return data_dir, params, steps, usable


def candidate_pool(data_dir):
    """The primers of `step3_df.csv`, in file order."""
    with open(os.path.join(data_dir, "step3_df.csv"), newline="") as handle:
        return [row["primer"] for row in csv.DictReader(handle)]


def csv_rows(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


OPTIMIZE_OUTPUT_PREFIX = "step4_"
FAILURE_FILE = "design_failure.json"


def _optimize_outputs(data_dir):
    return sorted(
        name
        for name in os.listdir(data_dir)
        if name.startswith(OPTIMIZE_OUTPUT_PREFIX) or name == FAILURE_FILE
    )


def run_optimize(design, data_dir, params, size, state, limits, retry):
    """`optimize` at one requested size, with its outputs moved to `sets/nNN/`.

    A design that fails is a result. Its `design_failure.json` is kept beside
    whatever else the step wrote, and the record says so.
    """
    from neoswga.core.delivered_set import read_delivered_set

    name = f"optimize_n{size:02d}"
    set_dir = os.path.join(data_dir, "sets", f"n{size:02d}")
    results_csv = os.path.join(set_dir, "step4_improved_df.csv")
    failure_json = os.path.join(set_dir, FAILURE_FILE)

    # Outputs left in the data directory by a run that was interrupted before
    # they were moved belong to no size; set them aside instead of reading them.
    leftovers = _optimize_outputs(data_dir)
    if leftovers:
        aside = os.path.join(data_dir, "sets", f"interrupted_{int(time.time())}")
        os.makedirs(aside, exist_ok=True)
        for leftover in leftovers:
            shutil.move(os.path.join(data_dir, leftover), os.path.join(aside, leftover))

    # `optimize` runs at the same worker count for every design, whatever the
    # pool steps used, so no delivered set differs from another by that.
    sized = dict(params, num_primers=size, target_set_size=size, cpus=DESIGN_PARAMS["cpus"])
    params_path = os.path.join(data_dir, f"params_n{size:02d}.json")
    with open(params_path, "w") as handle:
        json.dump(sized, handle, indent=2)
        handle.write("\n")

    def complete():
        return (
            os.path.exists(results_csv)
            or os.path.exists(failure_json)
            or bool(_optimize_outputs(data_dir))
        )

    record = run_step(
        name,
        cli_command("optimize", params_path, "--seed", str(OPTIMIZE_SEED)),
        data_dir,
        state,
        limits,
        retry,
        complete,
        design["id"],
        cpus=sized["cpus"],
    )

    produced = _optimize_outputs(data_dir)
    if produced:
        os.makedirs(set_dir, exist_ok=True)
        for produced_name in produced:
            shutil.move(os.path.join(data_dir, produced_name), os.path.join(set_dir, produced_name))
        manifest = os.path.join(data_dir, "run_manifest.json")
        if os.path.exists(manifest):
            shutil.copy2(manifest, os.path.join(set_dir, "run_manifest.json"))

    result = {
        "requested_size": size,
        "seed": OPTIMIZE_SEED,
        "step": record,
        "set_dir": os.path.relpath(set_dir, REPO_ROOT),
        "design_failure": None,
        "delivered_size": None,
        "delivered_primers": None,
        "alternative_sets": None,
        "optimizer_summary": None,
    }
    if os.path.exists(failure_json):
        with open(failure_json) as handle:
            result["design_failure"] = json.load(handle)
    if os.path.exists(results_csv):
        delivered = read_delivered_set(results_csv, 0)
        result["delivered_primers"] = list(delivered.primers)
        result["delivered_size"] = len(delivered.primers)
        result["alternative_sets"] = list(delivered.available)
        summary = os.path.join(set_dir, "step4_improved_df_summary.json")
        if os.path.exists(summary):
            with open(summary) as handle:
                result["optimizer_summary"] = json.load(handle)
    if result["delivered_primers"]:
        result["outcome"] = "delivered"
    elif result["design_failure"] is not None:
        result["outcome"] = "design_failed"
    else:
        result["outcome"] = f"no_set ({record['status']})"
    return result


# ----------------------------------------------------------------------
# Evaluation: the Phase 2 function, in a watched child
# ----------------------------------------------------------------------


def eval_prefix(workdir, key):
    return os.path.join(workdir, "eval", key)


def build_eval_tables(panel, workdir, limits, retry):
    """K-mer tables for every reference under `eval/`, by one `count-kmers`.

    These prefixes hold a table and no position index. A reference outside a
    design is read through them: positions by scanning, the weighted load from
    the table.
    """
    data_dir = os.path.join(workdir, "eval")
    os.makedirs(data_dir, exist_ok=True)
    design = {"id": "eval_tables", "fg": panel["targets"], "bg": panel["hosts"]}
    params = design_params(design, panel, data_dir)
    params.pop("genome_gc", None)
    digest = params_digest(params)
    state = load_state(data_dir)
    if state.get("params_digest") not in (None, digest):
        raise ValueError(f"{data_dir} was built under different parameters; it is not overwritten.")
    state["params_digest"] = digest
    params_path = os.path.join(data_dir, "params.json")
    with open(params_path, "w") as handle:
        json.dump(params, handle, indent=2)
        handle.write("\n")
    save_state(data_dir, state)
    return run_step(
        "count-kmers",
        cli_command("count-kmers", params_path),
        data_dir,
        state,
        limits,
        retry,
        lambda: pool_artifacts_present("count-kmers", params, data_dir),
        "eval_tables",
        cpus=params["cpus"],
    )


def build_counts_only_table(panel, key, workdir, limits, retry, cpus):
    """The k-mer table of one counts-only reference, by the pipeline's own step.

    Built the same way as the other seven, through `count-kmers`, so the table
    is the same kind of artifact read by the same door. Two things are
    deliberate.

    It has its own data directory and its own state, so adding a reference here
    does not change the parameter digest of `eval/` and the tables already
    there are neither refused nor recounted. The table itself still lands on
    the `eval/` prefix, which is where the evaluation looks for it.

    The reference is named as the FOREGROUND of that step with no background.
    `count-kmers` counts what it is given and reads no role, and the foreground
    list is the one that cannot be empty; calling it a background instead would
    have meant naming a foreground as well and recounting a table that is
    already correct.
    """
    genome = panel["genomes"][key]
    if key not in set(panel.get("counts_only_hosts", ())):
        raise ValueError(f"{key!r} is not a counts-only reference of this panel")
    data_dir = os.path.join(workdir, "eval_counts", key)
    os.makedirs(data_dir, exist_ok=True)
    prefix = eval_prefix(workdir, key)
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    params = dict(DESIGN_PARAMS)
    params["cpus"] = int(cpus)
    params["data_dir"] = data_dir
    params["fg_genomes"] = [genome["path"]]
    params["fg_prefixes"] = [prefix]
    params["fg_seq_lengths"] = [genome["length"]]
    params["bg_genomes"] = []
    params["bg_prefixes"] = []
    params["bg_seq_lengths"] = []
    params["fg_circular"] = False
    params.pop("genome_gc", None)
    digest = params_digest(params)
    state = load_state(data_dir)
    if state.get("params_digest") not in (None, digest):
        raise ValueError(f"{data_dir} was built under different parameters; it is not overwritten.")
    state["params_digest"] = digest
    params_path = os.path.join(data_dir, "params.json")
    with open(params_path, "w") as handle:
        json.dump(params, handle, indent=2)
        handle.write("\n")
    save_state(data_dir, state)
    record = run_step(
        "count-kmers",
        cli_command("count-kmers", params_path),
        data_dir,
        state,
        limits,
        retry,
        lambda: pool_artifacts_present("count-kmers", params, data_dir),
        f"counts_only_table {key}",
        cpus=params["cpus"],
    )
    from neoswga.core import kmer_tables

    files = kmer_tables.table_files(prefix, K) if kmer_tables.table_exists(prefix, K) else []
    return dict(
        record,
        key=key,
        prefix=os.path.relpath(prefix, REPO_ROOT),
        genome=os.path.relpath(genome["path"], REPO_ROOT),
        length_bp=genome["length"],
        records=genome.get("records"),
        table_files=[
            {
                "path": os.path.relpath(path, REPO_ROOT),
                "bytes": os.path.getsize(path),
            }
            for path in files
        ],
        table_bytes=sum(os.path.getsize(path) for path in files),
    )


def _position_index_exists(prefix):
    """Whether a position index file exists for this prefix at any k.

    The same question `reference_panel_evaluation._index_exists` asks, asked
    the same way: a directory listing, because the file itself is only ever
    read through `position_index`.
    """
    directory = os.path.dirname(prefix) or "."
    stem = os.path.basename(prefix)
    try:
        names = os.listdir(directory)
    except OSError:
        return False
    return any(name.startswith(stem + "_") and name.endswith("mer_positions.h5") for name in names)


def reference_plan(design, data_dir, panel, workdir, offdesign_host_scan, keys=None):
    """One entry per reference: which prefix answers for it, and by what route.

    Prefix and genome travel together on each entry. A reference in the design
    uses the design's prefix; any other uses the `eval/` prefix and a scan.
    `offdesign_host_scan=False` sends a host outside the design down the
    counts route, where sites are exact and coverage is unavailable.

    A counts-only reference takes the counts route whatever that flag says, and
    is refused outright if it somehow reached a design. `keys` narrows the plan
    to those references; it does not change any reference's route.
    """
    in_design = set(design["fg"]) | set(design["bg"])
    counts_only = set(panel.get("counts_only_hosts", ()))
    wanted = panel["targets"] + panel["hosts"] + list(panel.get("counts_only_hosts", ()))
    if keys is not None:
        unknown = [key for key in keys if key not in wanted]
        if unknown:
            raise ValueError(f"{unknown} are not references of this panel")
        wanted = [key for key in wanted if key in set(keys)]
    plan = []
    for key in wanted:
        genome = panel["genomes"][key]
        role = genome["role"]
        inside = key in in_design
        scan = True
        if key in counts_only:
            if inside:
                raise ValueError(
                    f"{key!r} is a counts-only reference and appears in design "
                    f"{design['id']!r}; that design should not exist"
                )
            # Not `and not offdesign_host_scan`: no flag reaches this.
            scan = False
            # `scan=False` is necessary and not sufficient. The Phase 2
            # function takes positions whenever an index exists for the
            # prefix, whatever `scan` says, so a counts-only reference is
            # only really on the counts route while its prefix has no index.
            # Nothing here builds one; this is the check that nothing did.
            if _position_index_exists(eval_prefix(workdir, key)):
                raise ValueError(
                    f"{eval_prefix(workdir, key)!r} has a position index. A "
                    f"counts-only reference must have none: the evaluation would "
                    f"read positions from it whatever `scan` says, which is the "
                    f"route this reference is too large for."
                )
        elif not inside and role == "host" and not offdesign_host_scan:
            scan = False
        if inside:
            route = "design index, scan for anything absent"
        elif scan:
            route = "scan of the FASTA"
        else:
            route = "k-mer counts; no positions"
        plan.append(
            {
                "key": key,
                "prefix": os.path.join(data_dir, key) if inside else eval_prefix(workdir, key),
                "genome": genome["path"],
                "length": genome["length"],
                "role": role,
                "circular": role == "target",
                "scan": scan,
                "in_design": inside,
                "held_out": key == design.get("held_out"),
                "route": route,
            }
        )
    return plan


def conditions_mapping():
    """The reaction every design here ran under, as the fields a builder reads."""
    return {
        key: DESIGN_PARAMS[key]
        for key in (
            "polymerase",
            "reaction_temp",
            "dmso_percent",
            "betaine_m",
            "trehalose_m",
            "na_conc",
            "mg_conc",
        )
    }


def evaluate_panels(panels, references, with_weighted_load, out_path, log_path, limits, label):
    """Evaluate labelled panels in one watched child; return its record and output."""
    spec = {
        "panels": panels,
        "references": references,
        "conditions": conditions_mapping() if with_weighted_load else None,
    }
    spec_path = out_path + ".spec.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(spec_path, "w") as handle:
        json.dump(spec, handle, indent=2)
        handle.write("\n")
    if os.path.exists(out_path):
        with open(out_path) as handle:
            cached = json.load(handle)
        if cached.get("spec") == spec:
            # The cost of the run that produced it is kept beside the output.
            record = {"status": STATUS_OK}
            if os.path.exists(out_path + ".step.json"):
                with open(out_path + ".step.json") as handle:
                    record = json.load(handle)
            return dict(record, reused=True), cached
    print(f"  {label}: evaluating {len(panels)} panel(s)", flush=True)
    record = run_watched(
        [sys.executable, os.path.abspath(__file__), "--worker", spec_path, out_path],
        log_path,
        limits,
        label,
    )
    record["reused"] = False
    with open(out_path + ".step.json", "w") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    check_halt(record, limits, f"{label} evaluation")
    if record["status"] != STATUS_OK or not os.path.exists(out_path):
        return record, None
    with open(out_path) as handle:
        return record, json.load(handle)


def worker_main(spec_path, out_path):
    """Child entry: call the Phase 2 function once per panel and write the result."""
    from types import SimpleNamespace

    from neoswga.core.coverage import polymerase_extension_reach
    from neoswga.core.reaction_conditions import build_reaction_conditions
    from neoswga.core.reference_panel_evaluation import ReferenceSpec, evaluate_reference_panel

    with open(spec_path) as handle:
        spec = json.load(handle)
    conditions = None
    if spec["conditions"] is not None:
        # Built from this mapping alone, so nothing a previous command left in
        # the `parameter` module can reach the chemistry.
        conditions = build_reaction_conditions(
            SimpleNamespace(**spec["conditions"]), from_mapping_only=True
        )
    reach = polymerase_extension_reach(DESIGN_PARAMS["polymerase"], coverage_metric="realistic")
    specs = [
        ReferenceSpec(
            prefix=ref["prefix"],
            genome=ref["genome"],
            length=int(ref["length"]),
            role=ref["role"],
            circular=bool(ref["circular"]),
            scan=bool(ref["scan"]),
        )
        for ref in spec["references"]
    ]
    out = {"spec": spec, "extension_reach_bp": int(reach), "panels": []}
    for panel in spec["panels"]:
        started = time.time()
        assessment = evaluate_reference_panel(
            panel["primers"], specs, conditions, extension=int(reach)
        )
        out["panels"].append(
            {
                "label": panel["label"],
                "seconds": round(time.time() - started, 2),
                "assessment": assessment.as_dict(),
            }
        )
    with open(out_path + ".tmp", "w") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")
    os.replace(out_path + ".tmp", out_path)
    return 0


def by_key(assessment, references):
    """Re-key one assessment from prefixes to manifest keys.

    The Phase 2 function keys its blocks by prefix. A reader comparing designs
    wants the genome, so the prefix is kept inside each record and the manifest
    key becomes the dictionary key. Nothing is recomputed.
    """
    names = {ref["prefix"]: ref for ref in references}

    def rekey(block):
        out = {}
        for prefix, record in block.items():
            ref = names[prefix]
            out[ref["key"]] = dict(
                record,
                prefix=os.path.relpath(prefix, REPO_ROOT),
                genome=os.path.relpath(record["genome"], REPO_ROOT) if record["genome"] else None,
                in_design=ref["in_design"],
                held_out=ref["held_out"],
                route=ref["route"],
            )
        return out

    pairs = []
    for pair in assessment["target_host_pairs"]:
        target = assessment["per_target"][pair["target"]]
        host = assessment["per_host"][pair["host"]]
        pairs.append(
            dict(
                pair,
                target=names[pair["target"]]["key"],
                host=names[pair["host"]]["key"],
                weighted_selectivity_density=weighted_density(target, host),
            )
        )
    return {
        "primers": assessment["primers"],
        "extension_reach_bp": assessment["extension_reach_bp"],
        "per_target": rekey(assessment["per_target"]),
        "per_host": rekey(assessment["per_host"]),
        "target_host_pairs": pairs,
        "worst_target_coverage": assessment["worst_target_coverage"],
        "worst_host_selectivity_density": assessment["worst_host_selectivity_density"],
        "notes": assessment["notes"],
    }


def weighted_density(target, host):
    """Target against host density from the two occupancy-weighted loads.

    DERIVED HERE, not returned by the Phase 2 function, whose pair density is
    formed from exact site counts. This one is formed from the weighted loads
    that function did return, with `selectivity_density_from_loads`, and is the
    quantity the optimizer reports as `selectivity_density` in occupancy mode.
    It is unavailable when either load is.
    """
    from neoswga.core.selectivity import selectivity_density_from_loads

    for record, role in ((target, "target"), (host, "host")):
        load = record["weighted_site_load"]
        if load["value"] is None:
            return {
                "value": None,
                "units": "ratio",
                "basis": None,
                "unavailable": f"the {role} has no weighted load: {load['unavailable']}",
            }
    fg = float(target["weighted_site_load"]["value"])
    bg = float(host["weighted_site_load"]["value"])
    return {
        "value": float(
            selectivity_density_from_loads(
                fg, float(target["length_bp"]), bg, float(host["length_bp"])
            )
        ),
        "units": "ratio",
        "basis": (
            f"derived by this script: weighted load {fg:.3f} per {target['length_bp']} bp "
            f"against {bg:.3f} per {host['length_bp']} bp"
        ),
        "unavailable": None,
    }


def headline(evaluation):
    """The per-reference and per-pair figures of one evaluation, without the
    per-primer maps. Each value is the `Measurement` as the function returned it.
    """
    fields = ("sites", "sites_per_mb", "coverage", "max_gap_bp", "weighted_site_load")
    return {
        "per_target": {
            key: {name: record[name] for name in fields}
            for key, record in evaluation["per_target"].items()
        },
        "per_host": {
            key: {name: record[name] for name in fields}
            for key, record in evaluation["per_host"].items()
        },
        "target_host_pairs": [
            {
                "target": pair["target"],
                "host": pair["host"],
                "selectivity_density": pair["selectivity_density"],
                "weighted_selectivity_density": pair["weighted_selectivity_density"],
                "zero_host_sites": pair["zero_host_sites"],
            }
            for pair in evaluation["target_host_pairs"]
        ],
        "worst_target_coverage": evaluation["worst_target_coverage"],
        "worst_host_selectivity_density": evaluation["worst_host_selectivity_density"],
    }


def spread(values, what):
    """Minimum, median and maximum over seeds, or unavailable with the reason.

    One unavailable member makes the whole summary unavailable: a spread over
    the seeds that happened to answer is not the spread.
    """
    missing = sum(1 for value in values if value is None)
    if missing or not values:
        return {
            "n": len(values),
            "min": None,
            "median": None,
            "max": None,
            "unavailable": f"{missing} of {len(values)} seed(s) have no measured {what}",
        }
    return {
        "n": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "unavailable": None,
    }


def summarise_seeds(evaluations):
    """Per-strain coverage and per-pair density across the C2 seeds."""
    if not evaluations:
        return None
    first = evaluations[0]
    coverage = {
        key: spread([e["per_target"][key]["coverage"]["value"] for e in evaluations], "coverage")
        for key in first["per_target"]
    }
    host_sites = {
        key: spread([e["per_host"][key]["sites_per_mb"]["value"] for e in evaluations], "sites")
        for key in first["per_host"]
    }
    density = {}
    weighted = {}
    for index, pair in enumerate(first["target_host_pairs"]):
        name = f"{pair['target']}|{pair['host']}"
        density[name] = spread(
            [e["target_host_pairs"][index]["selectivity_density"]["value"] for e in evaluations],
            "selectivity density",
        )
        # A pair with no host sites carries the package's ceiling value, not a
        # measured ratio; the count says how many seeds that applies to.
        density[name]["seeds_with_zero_host_sites"] = sum(
            1 for e in evaluations if e["target_host_pairs"][index]["zero_host_sites"]
        )
        weighted[name] = spread(
            [
                e["target_host_pairs"][index]["weighted_selectivity_density"]["value"]
                for e in evaluations
            ],
            "weighted selectivity density",
        )
    return {
        "coverage_by_target": coverage,
        "host_sites_per_mb": host_sites,
        "selectivity_density_by_pair": density,
        "weighted_selectivity_density_by_pair": weighted,
    }


def random_panels(pool, size, seeds):
    """One uniform draw of `size` primers per seed, from the sorted pool."""
    ordered = sorted(pool)
    panels = []
    for seed in range(seeds):
        rng = random.Random(seed)
        panels.append({"label": f"seed_{seed:02d}", "primers": rng.sample(ordered, size)})
    return panels


# ----------------------------------------------------------------------
# The run
# ----------------------------------------------------------------------


def compare_with_example_pool(data_dir):
    """Whether this design's candidate pool is the prepared example's pool.

    Read as CSV only. The recorded reference figure was measured on that pool,
    so a baseline is the same pool only if the candidate sets are equal.
    """
    example = os.path.join(REPO_ROOT, "examples", "wolbachia_pool_design", "work", "step3_df.csv")
    if not os.path.exists(example):
        return {"compared": False, "reason": f"{os.path.relpath(example, REPO_ROOT)} is absent"}
    with open(example, newline="") as handle:
        theirs = [row["primer"] for row in csv.DictReader(handle)]
    ours = candidate_pool(data_dir)
    return {
        "compared": True,
        "example_pool": os.path.relpath(example, REPO_ROOT),
        "example_candidates": len(theirs),
        "candidates": len(ours),
        "shared": len(set(ours) & set(theirs)),
        "same_set": set(ours) == set(theirs),
        "same_order": ours == theirs,
    }


def file_digest(path):
    return sha256_of(path) if os.path.exists(path) else None


# ----------------------------------------------------------------------
# The counts route: one agreement check, then the counts-only references
# ----------------------------------------------------------------------


def load_results(output):
    """An existing results file, refused if it holds no design.

    A pass that only adds figures must not be the thing that creates the file:
    an empty or absent results file means the designs are not there to measure,
    and writing one here would look like a run that produced them.
    """
    if not os.path.exists(output):
        raise ValueError(
            f"{os.path.relpath(output, REPO_ROOT)} does not exist. This pass adds "
            f"figures to the designs an earlier run recorded; it builds none."
        )
    with open(output) as handle:
        results = json.load(handle)
    if not results.get("designs"):
        raise ValueError(f"{os.path.relpath(output, REPO_ROOT)} records no design.")
    return results


def _same_panels(raw_path, panels):
    """Whether these panels are the ones that output was produced from.

    True or False when the earlier output can be read, None when it cannot: an
    absent file says the draws were not checked, which is not the same as their
    being wrong.
    """
    if not os.path.exists(raw_path):
        return None
    with open(raw_path) as handle:
        spec = json.load(handle).get("spec") or {}
    was = spec.get("panels")
    if not was:
        return None
    return [(p["label"], list(p["primers"])) for p in was] == [
        (p["label"], list(p["primers"])) for p in panels
    ]


def design_of(entry):
    """The design dictionary an earlier run planned, back from its record."""
    return {
        "id": entry["id"],
        "group": entry["group"],
        "fg": entry["fg"],
        "bg": entry["bg"],
        "held_out": entry["held_out"],
    }


def counts_route_check(args, panel, workdir, limits, output):
    """Measure one host both ways on one recorded panel and compare the sites.

    The counts route and the position route answer the same question about
    exact binding sites, at very different cost, and only one of them can also
    answer a positional one. This checks that claim on this data before any
    figure is read through the counts route alone: the host is measured again
    by counts, and the site count is compared with the one already recorded for
    the same panel through positions.

    A disagreement is the finding. Nothing is reconciled here.
    """
    results = load_results(output)
    host = args.counts_route_check_host
    entries = [e for e in results["designs"] if e["id"] == args.counts_route_check]
    if not entries:
        raise ValueError(
            f"no design {args.counts_route_check!r} in "
            f"{os.path.relpath(output, REPO_ROOT)}; it records "
            + ", ".join(e["id"] for e in results["designs"])
        )
    entry = entries[0]
    sized = [s for s in entry["sizes"] if s["requested_size"] == args.counts_route_check_size]
    if not sized or not sized[0].get("evaluation"):
        raise ValueError(
            f"{entry['id']} has no evaluated set at requested size "
            f"{args.counts_route_check_size}"
        )
    sized = sized[0]
    recorded = sized["evaluation"]["per_host"].get(host)
    if recorded is None:
        raise ValueError(f"{entry['id']} was not evaluated against {host!r}")

    design = design_of(entry)
    data_dir = os.path.join(REPO_ROOT, entry["data_dir"])
    # The same host, the same genome and the same panel, measured by counts.
    #
    # `scan=False` alone does not force that route: the Phase 2 function takes
    # positions whenever an index exists for the prefix, whatever `scan` says,
    # and the design's own prefix has one for every host it was designed
    # against. So the counts route is reached by naming the `eval/` prefix,
    # which holds a table and no index. That is the same prefix a counts-only
    # reference is read through, which is the route being checked.
    counted = [
        dict(
            ref,
            prefix=eval_prefix(workdir, host),
            scan=False,
            route="k-mer counts; no positions",
        )
        for ref in reference_plan(design, data_dir, panel, workdir, True, keys=[host])
    ]
    target = reference_plan(design, data_dir, panel, workdir, True, keys=[panel["reference"]])
    plan = target + counted
    references = [{k: v for k, v in ref.items() if k != "route"} for ref in plan]
    out_path = os.path.join(
        workdir, "counts_route_check", f"{entry['id']}_n{sized['requested_size']:02d}.json"
    )
    record, out = evaluate_panels(
        [{"label": "set_0", "primers": sized["delivered_primers"]}],
        references,
        False,
        out_path,
        os.path.join(workdir, "counts_route_check", "logs", f"{entry['id']}.log"),
        limits,
        f"counts-route check {entry['id']} n={sized['requested_size']}",
    )
    check = {
        "design": entry["id"],
        "requested_size": sized["requested_size"],
        "delivered_size": sized["delivered_size"],
        "host": host,
        "panel": sized["delivered_primers"],
        "step": record,
        "positions_route": {
            "route": recorded["route"],
            "site_source": recorded["site_source"],
            "sites": recorded["sites"],
            "coverage": recorded["coverage"],
        },
        "counts_route": None,
        "agree": None,
    }
    if out is not None:
        measured = by_key(out["panels"][0]["assessment"], plan)["per_host"][host]
        check["counts_route"] = {
            "route": "k-mer counts; no positions",
            "site_source": measured["site_source"],
            "sites": measured["sites"],
            "coverage": measured["coverage"],
        }
        left = recorded["sites"]["value"]
        right = measured["sites"]["value"]
        check["agree"] = (
            None if left is None or right is None else bool(float(left) == float(right))
        )
        check["difference"] = None if left is None or right is None else float(right) - float(left)
    results.setdefault("counts_route_checks", []).append(check)
    with open(output + ".tmp", "w") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")
    os.replace(output + ".tmp", output)
    print(json.dumps(check, indent=2))
    return 0 if check["agree"] else 1


def run_counts_only_pass(args, panel, workdir, limits, output):
    """Measure every recorded panel against the counts-only references.

    Reads the results file, adds figures to it, and writes it back. Nothing
    else in it is touched: the designs, their step records and the evaluations
    an earlier run wrote are carried through as they were read.

    Only the targets and the counts-only references are in each plan. The hosts
    the design was evaluated against already have their figures, and measuring
    them again would scan a 144 Mb reference once per panel for an answer the
    file already holds.
    """
    results = load_results(output)
    keys = list(panel["counts_only_hosts"])
    if not keys:
        raise ValueError(
            "no counts-only reference was named. Pass --counts-only-hosts KEY, or "
            f"name a host longer than {COUNTS_ONLY_ABOVE_BP:,} bp."
        )
    designs_before = [e["id"] for e in results["designs"]]
    steps_before = {e["id"]: sorted(e.get("steps", {})) for e in results["designs"]}
    # The reader of the tables asks the top-level key, which only a full run
    # sets. Without it a merged-in reference is in the file and absent from
    # every table, which reads as a reference that was not measured.
    results["counts_only_hosts"] = keys
    # And into the one mapping a reader asks for any reference's length, role
    # and record count. Leaving a counts-only reference out of it while it
    # appears in `counts_only_hosts` is a reference named in one place and
    # described in none.
    for key in keys:
        genome = panel["genomes"][key]
        results["genomes"][key] = {k: v for k, v in genome.items() if k != "path"} | {
            "path": os.path.relpath(genome["path"], REPO_ROOT)
        }
    # A halt recorded by an earlier attempt at this pass is not a result of
    # this one. It is cleared here and re-set below only if this pass halts.
    results.pop("counts_only_halted", None)

    def write():
        results["written"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        results["free_percent_at_write"] = free_percent()
        after = [e["id"] for e in results["designs"]]
        if after != designs_before:
            raise RuntimeError(f"the designs changed: {designs_before} became {after}")
        for entry in results["designs"]:
            if sorted(entry.get("steps", {})) != steps_before[entry["id"]]:
                raise RuntimeError(f"the step records of {entry['id']} changed")
        with open(output + ".tmp", "w") as handle:
            json.dump(results, handle, indent=2)
            handle.write("\n")
        os.replace(output + ".tmp", output)

    block = results.setdefault(
        "counts_only",
        {
            "note": (
                "References measured by k-mer counts alone. Sites and densities are "
                "exact; coverage, the gap figures and host coverage are unavailable, "
                "with the reason carried on each figure. No design has one of these "
                "as a background: they were never counted into a design, never "
                "filtered against and never scanned for positions."
            ),
            "above_bp": COUNTS_ONLY_ABOVE_BP,
            "hosts": {},
            "tables": {},
        },
    )
    block["hosts"] = {
        key: {k: v for k, v in panel["genomes"][key].items() if k != "path"}
        | {"path": os.path.relpath(panel["genomes"][key]["path"], REPO_ROOT)}
        for key in keys
    }
    write()

    try:
        for key in keys:
            # One at a time. A 3.3 Gb count is the largest single step of this
            # stage and nothing else runs beside it.
            block["tables"][key] = build_counts_only_table(
                panel, key, workdir, limits, args.retry, args.counts_only_cpus
            )
            write()
            if block["tables"][key]["status"] != STATUS_OK:
                raise RunHalted(
                    f"the {key} table could not be counted; see "
                    f"{block['tables'][key].get('log')}"
                )
        if args.counts_only_tables_only:
            write()
            for key in keys:
                record = block["tables"][key]
                peak = record.get("peak_tree_rss_bytes")
                print(
                    f"{key}: {record['status']} in {record['wall_seconds']} s, peak "
                    + ("not sampled" if peak is None else f"{peak / 2**30:.2f} GiB")
                    + f", table {record['table_bytes'] / 2**20:.1f} MiB, free "
                    f"{record['free_percent_before']} -> {record['free_percent_lowest']} -> "
                    f"{record['free_percent_after']} percent"
                )
            return 0

        for entry in results["designs"]:
            design = design_of(entry)
            data_dir = os.path.join(REPO_ROOT, entry["data_dir"])
            plan = reference_plan(
                design, data_dir, panel, workdir, False, keys=panel["targets"] + keys
            )
            references = [{k: v for k, v in ref.items() if k != "route"} for ref in plan]
            entry["counts_only_references"] = [
                {k: ref[k] for k in ("key", "role", "in_design", "held_out", "route")}
                for ref in plan
            ]
            for sized in entry["sizes"]:
                if sized.get("outcome") != "delivered":
                    sized["counts_only"] = None
                    sized["counts_only_unavailable"] = f"no delivered set: {sized.get('outcome')}"
                    write()
                    continue
                set_dir = os.path.join(REPO_ROOT, sized["set_dir"])
                label = f"{entry['id']} n={sized['requested_size']} counts-only"
                record, out = evaluate_panels(
                    [{"label": "set_0", "primers": sized["delivered_primers"]}],
                    references,
                    True,
                    os.path.join(set_dir, "counts_only_evaluation.json"),
                    os.path.join(
                        data_dir, "logs", f"counts_only_n{sized['requested_size']:02d}.log"
                    ),
                    limits,
                    label,
                )
                sized["counts_only_step"] = record
                if out is None:
                    sized["counts_only"] = None
                    sized["counts_only_unavailable"] = (
                        f"the evaluation child ended with status {record['status']}; "
                        f"see {record.get('log')}"
                    )
                else:
                    sized["counts_only"] = by_key(out["panels"][0]["assessment"], plan)
                write()

                recorded = sized.get("c2")
                if not recorded or not recorded.get("panels"):
                    sized["c2_counts_only"] = None
                    sized["c2_counts_only_unavailable"] = (
                        "no C2 panels were recorded for this design and size"
                    )
                    write()
                    continue
                pool = candidate_pool(data_dir)
                panels = random_panels(pool, recorded["panel_size"], recorded["seeds"])
                # The same draws, or the two blocks describe different panels
                # and reading them side by side is wrong. `random_panels` is
                # deterministic in the pool, the size and the seed, so this
                # holds by construction; it is checked against the panels the
                # earlier pass actually evaluated, which its own output carries.
                c2_same = _same_panels(os.path.join(REPO_ROOT, recorded["raw"]), panels)
                if c2_same is False:
                    raise RunHalted(
                        f"the C2 draws for {entry['id']} at size {recorded['panel_size']} "
                        f"do not match the ones already recorded in {recorded['raw']}; "
                        f"the candidate pool must have changed"
                    )
                record, out = evaluate_panels(
                    panels,
                    references,
                    True,
                    os.path.join(set_dir, "c2_counts_only_evaluation.json"),
                    os.path.join(
                        data_dir, "logs", f"c2_counts_only_n{sized['requested_size']:02d}.log"
                    ),
                    limits,
                    f"{label} C2 at delivered size {recorded['panel_size']}",
                )
                c2 = {
                    "panel_size": recorded["panel_size"],
                    "seeds": recorded["seeds"],
                    "pool_candidates": len(pool),
                    "same_panels_as_recorded_c2": c2_same,
                    "step": record,
                    "raw": os.path.relpath(
                        os.path.join(set_dir, "c2_counts_only_evaluation.json"), REPO_ROOT
                    ),
                }
                if out is None:
                    c2["panels"] = None
                    c2["summary"] = None
                    c2["unavailable"] = f"the evaluation child ended with status {record['status']}"
                else:
                    evaluations = [
                        by_key(panel_out["assessment"], plan) for panel_out in out["panels"]
                    ]
                    c2["panels"] = [
                        dict(headline(evaluation), label=panel_out["label"])
                        for evaluation, panel_out in zip(evaluations, out["panels"], strict=True)
                    ]
                    c2["summary"] = summarise_seeds(evaluations)
                sized["c2_counts_only"] = c2
                write()
    except RunHalted as halted:
        results["counts_only_halted"] = str(halted)
        write()
        print(f"HALTED: {halted}", flush=True)
        return EXIT_HALTED

    write()
    print(f"Merged the counts-only figures into {os.path.relpath(output, REPO_ROOT)}")
    return 0


def cpus_check(design, panel, workdir, limits, retry, cpus):
    """Rebuild one candidate pool at another worker count and compare it.

    The pool is rebuilt in a separate data directory, and `step2_df.csv` and
    `step3_df.csv` are compared byte for byte with the ones in the design's
    own directory. A file absent on either side is `compared: False` with the
    reason, never a match.
    """
    baseline_dir = os.path.join(workdir, "pools", design["id"])
    check_dir = os.path.join(workdir, "cpus_check", f"{design['id']}__cpus{cpus}")
    _dir, _params, steps, usable = build_pool(
        design, panel, workdir, limits, retry, cpus, data_dir=check_dir
    )
    baseline_state = load_state(baseline_dir)["steps"] if os.path.isdir(baseline_dir) else {}
    files = {}
    for name in ("step2_df.csv", "step3_df.csv"):
        ours = file_digest(os.path.join(check_dir, name))
        theirs = file_digest(os.path.join(baseline_dir, name))
        if ours is None or theirs is None:
            files[name] = {
                "compared": False,
                "reason": "absent from "
                + " and ".join(
                    label
                    for label, digest in (("the check", ours), ("the baseline", theirs))
                    if digest is None
                ),
            }
        else:
            files[name] = {
                "compared": True,
                "identical": ours == theirs,
                "sha256_check": ours,
                "sha256_baseline": theirs,
                "bytes": os.path.getsize(os.path.join(check_dir, name)),
            }
    return {
        "design": design["id"],
        "cpus_check": cpus,
        "check_dir": os.path.relpath(check_dir, REPO_ROOT),
        "baseline_dir": os.path.relpath(baseline_dir, REPO_ROOT),
        "pool_built": usable,
        "steps_check": steps,
        "steps_baseline": {
            name: baseline_state.get(name)
            for name in ("count-kmers", "filter", "prepare-candidates")
        },
        "files": files,
    }


def load_cpus_checks(workdir):
    path = os.path.join(workdir, "cpus_check.json")
    if not os.path.exists(path):
        return []
    with open(path) as handle:
        return json.load(handle)


def run_cpus_check(args, panel, workdir, limits, c1_host):
    matches = [d for d in plan_designs(panel, c1_host) if d["id"] == args.cpus_check]
    if not matches:
        raise ValueError(f"--cpus-check: no design has the id {args.cpus_check!r}")
    if not args.pooled_cpus:
        raise ValueError("--cpus-check needs --pooled-cpus, the worker count to compare")
    try:
        record = cpus_check(matches[0], panel, workdir, limits, args.retry, args.pooled_cpus)
    except RunHalted as halted:
        print(f"HALTED: {halted}", flush=True)
        return EXIT_HALTED
    checks = [
        c
        for c in load_cpus_checks(workdir)
        if (c["design"], c["cpus_check"]) != (record["design"], record["cpus_check"])
    ]
    checks.append(record)
    path = os.path.join(workdir, "cpus_check.json")
    with open(path + ".tmp", "w") as handle:
        json.dump(checks, handle, indent=2)
        handle.write("\n")
    os.replace(path + ".tmp", path)
    for name, outcome in record["files"].items():
        if outcome["compared"]:
            verdict = "identical" if outcome["identical"] else "DIFFERENT"
        else:
            verdict = "not compared: " + outcome["reason"]
        print(f"  {name}: {verdict}")
    print(f"Wrote {os.path.relpath(path, REPO_ROOT)}")
    return 0


def run(args):
    workdir = os.path.abspath(args.workdir)
    ignored = subprocess.run(
        ["git", "check-ignore", "-q", os.path.join(workdir, "probe")], cwd=REPO_ROOT, check=False
    )
    if ignored.returncode != 0:
        raise ValueError(f"{workdir} is not gitignored; generated tables must not be tracked.")
    os.makedirs(workdir, exist_ok=True)
    output = os.path.abspath(args.output or os.path.join(workdir, "results.json"))

    limits = {
        "rss_limit_bytes": int(args.rss_limit_gb * 2**30),
        "halt_peak_bytes": int(args.halt_peak_gb * 2**30),
        "stop_free_percent": args.stop_free_percent,
        "start_free_percent": args.start_free_percent,
        "start_wait_seconds": args.start_wait_seconds,
        "sample_seconds": args.sample_seconds,
        "min_free_disk_bytes": int(args.min_free_disk_gb * 2**30),
        "disk_path": workdir,
    }
    groups = [group.strip() for group in args.only.split(",")] if args.only else list(GROUPS)
    unknown = [group for group in groups if group not in GROUPS]
    if unknown:
        raise ValueError(f"unknown group(s) {unknown}; choose from {list(GROUPS)}")
    sizes = [int(size) for size in args.sizes.split(",")]
    hosts = [host.strip() for host in args.hosts.split(",") if host.strip()]
    counts_only = [key.strip() for key in (args.counts_only_hosts or "").split(",") if key.strip()]

    panel = load_panel(args.manifest, hosts, args.max_host_bp, counts_only)
    c1_host = args.c1_host or max(panel["hosts"], key=lambda key: panel["genomes"][key]["length"])
    if c1_host not in panel["hosts"]:
        raise ValueError(f"--c1-host {c1_host!r} is not among the hosts of this run")

    if args.cpus_check:
        return run_cpus_check(args, panel, workdir, limits, c1_host)
    if args.counts_route_check:
        return counts_route_check(args, panel, workdir, limits, output)
    if args.counts_only_pass:
        return run_counts_only_pass(args, panel, workdir, limits, output)
    if args.tables:
        print_tables(load_results(output))
        return 0

    results = {
        "script": os.path.relpath(os.path.abspath(__file__), REPO_ROOT),
        "manifest": os.path.relpath(os.path.abspath(args.manifest), REPO_ROOT),
        "workdir": os.path.relpath(workdir, REPO_ROOT),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": sys.version.split()[0],
        "git_sha": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True, check=False
        ).stdout.strip()
        or None,
        "k": K,
        "sizes_requested": sizes,
        "groups_requested": groups,
        "c2_seeds": args.seeds,
        "c2_note": (
            "Uniform draws from step3_df.csv at the delivered size. Not dimer-screened. "
            + (
                "Hosts outside the design were scanned."
                if args.c2_scan_offdesign_hosts
                else "Sites on a host outside the design come from k-mer counts, so "
                "coverage on that host is unavailable for these panels."
            )
        ),
        "design_parameters": DESIGN_PARAMS,
        "optimize_seed": OPTIMIZE_SEED,
        "limits": limits,
        "targets": panel["targets"],
        "reference": panel["reference"],
        "hosts": panel["hosts"],
        "counts_only_hosts": panel["counts_only_hosts"],
        "hosts_not_run": {
            key: "not named in --hosts for this run; nothing about it was measured"
            for key in panel["hosts_not_run"]
        },
        "genomes": {
            key: {k: v for k, v in genome.items() if k != "path"}
            | {"path": os.path.relpath(genome["path"], REPO_ROOT)}
            for key, genome in panel["genomes"].items()
        },
        "recorded_reference": RECORDED_REFERENCE,
        "eval_tables": None,
        "cpus_checks": load_cpus_checks(workdir),
        "designs": [],
        "halted": None,
    }

    def write():
        results["written"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        results["free_percent_at_write"] = free_percent()
        with open(output + ".tmp", "w") as handle:
            json.dump(results, handle, indent=2)
            handle.write("\n")
        os.replace(output + ".tmp", output)

    wanted = [
        d
        for d in plan_designs(panel, c1_host)
        if d["group"] in groups and (not args.match or args.match in d["id"])
    ]
    try:
        results["eval_tables"] = build_eval_tables(panel, workdir, limits, args.retry)
        write()
        if results["eval_tables"]["status"] != STATUS_OK:
            raise RunHalted("the evaluation tables could not be counted; see eval/logs")

        built = []
        for design in wanted:
            print(f"{design['id']}", flush=True)
            entry = {
                "id": design["id"],
                "group": design["group"],
                "label": (
                    "preliminary D3: the small hosts only, not the full host panel"
                    if design["group"] == "D3"
                    else design["group"]
                ),
                "fg": design["fg"],
                "bg": design["bg"],
                "held_out": design["held_out"],
                "sizes": [],
            }
            results["designs"].append(entry)
            pooled = len(design["fg"]) > 1
            pool_cpus = args.pooled_cpus if pooled and args.pooled_cpus else DESIGN_PARAMS["cpus"]
            data_dir, params, steps, usable = build_pool(
                design, panel, workdir, limits, args.retry, pool_cpus
            )
            entry["cpus"] = {
                "pool_steps_requested": pool_cpus,
                "pool_steps_ran_under": {name: rec.get("cpus") for name, rec in steps.items()},
                "optimize": DESIGN_PARAMS["cpus"],
            }
            entry["data_dir"] = os.path.relpath(data_dir, REPO_ROOT)
            entry["steps"] = steps
            entry["genome_gc_given"] = params.get("genome_gc")
            entry["step2_candidates"] = csv_rows(os.path.join(data_dir, "step2_df.csv"))
            entry["step3_candidates"] = csv_rows(os.path.join(data_dir, "step3_df.csv"))
            write()
            if not usable:
                entry["outcome"] = "the candidate pool was not built; see steps"
                continue
            if design["fg"] == [panel["reference"]] and design["bg"] == ["drosophila"]:
                entry["example_pool_comparison"] = compare_with_example_pool(data_dir)

            state = load_state(data_dir)
            references = reference_plan(design, data_dir, panel, workdir, True)
            entry["references"] = [
                {k: ref[k] for k in ("key", "role", "in_design", "held_out", "route")}
                for ref in references
            ]
            for size in sizes:
                sized = run_optimize(design, data_dir, params, size, state, limits, args.retry)
                entry["sizes"].append(sized)
                write()
                if sized["outcome"] != "delivered":
                    sized["evaluation"] = None
                    sized["evaluation_unavailable"] = f"no delivered set: {sized['outcome']}"
                    continue
                set_dir = os.path.join(REPO_ROOT, sized["set_dir"])
                record, out = evaluate_panels(
                    [{"label": "set_0", "primers": sized["delivered_primers"]}],
                    [{k: v for k, v in ref.items() if k != "route"} for ref in references],
                    True,
                    os.path.join(set_dir, "evaluation.json"),
                    os.path.join(data_dir, "logs", f"evaluate_n{size:02d}.log"),
                    limits,
                    f"{design['id']} n={size}",
                )
                sized["evaluation_step"] = record
                if out is None:
                    sized["evaluation"] = None
                    sized["evaluation_unavailable"] = (
                        f"the evaluation child ended with status {record['status']}; "
                        f"see {record.get('log')}"
                    )
                else:
                    sized["evaluation_seconds"] = out["panels"][0]["seconds"]
                    sized["evaluation"] = by_key(out["panels"][0]["assessment"], references)
                write()
            built.append((design, entry, data_dir))

        if "C2" in groups:
            for design, entry, data_dir in built:
                if not design["c2"]:
                    continue
                pool = candidate_pool(data_dir)
                references = reference_plan(
                    design, data_dir, panel, workdir, args.c2_scan_offdesign_hosts
                )
                for sized in entry["sizes"]:
                    if sized["outcome"] != "delivered":
                        sized["c2"] = None
                        sized["c2_unavailable"] = "no delivered set to match in size"
                        continue
                    size = sized["delivered_size"]
                    set_dir = os.path.join(REPO_ROOT, sized["set_dir"])
                    record, out = evaluate_panels(
                        random_panels(pool, size, args.seeds),
                        [{k: v for k, v in ref.items() if k != "route"} for ref in references],
                        True,
                        os.path.join(set_dir, "c2_evaluation.json"),
                        os.path.join(data_dir, "logs", f"c2_n{sized['requested_size']:02d}.log"),
                        limits,
                        f"{design['id']} C2 at delivered size {size}",
                    )
                    c2 = {
                        "panel_size": size,
                        "seeds": args.seeds,
                        "pool_candidates": len(pool),
                        "step": record,
                        "references": [
                            {k: ref[k] for k in ("key", "role", "in_design", "route")}
                            for ref in references
                        ],
                        "raw": os.path.relpath(
                            os.path.join(set_dir, "c2_evaluation.json"), REPO_ROOT
                        ),
                    }
                    if out is None:
                        c2["panels"] = None
                        c2["summary"] = None
                        c2["unavailable"] = (
                            f"the evaluation child ended with status {record['status']}"
                        )
                    else:
                        evaluations = [
                            by_key(panel_out["assessment"], references)
                            for panel_out in out["panels"]
                        ]
                        c2["panels"] = [
                            dict(headline(evaluation), label=panel_out["label"])
                            for evaluation, panel_out in zip(
                                evaluations, out["panels"], strict=True
                            )
                        ]
                        c2["summary"] = summarise_seeds(evaluations)
                    sized["c2"] = c2
                    write()
    except RunHalted as halted:
        results["halted"] = str(halted)
        write()
        print(f"HALTED: {halted}", flush=True)
        print(f"Wrote {os.path.relpath(output, REPO_ROOT)}")
        return EXIT_HALTED

    write()
    print_report(results)
    print(f"\nWrote {os.path.relpath(output, REPO_ROOT)}")
    return 0


# ----------------------------------------------------------------------
# The report
# ----------------------------------------------------------------------


def _value(measurement, fmt):
    if measurement is None or measurement.get("value") is None:
        return "n/a"
    return fmt.format(measurement["value"])


def print_report(results):
    targets = results["targets"]
    hosts = results["hosts"]
    print()
    print("Coverage per strain and selectivity density per host, set 0 of each design.")
    print("n/a is a figure that was not measured; the JSON carries the reason.")
    header = f"{'design':<44}{'req':>4}{'got':>4}"
    for key in targets:
        header += f"{key.replace('wolbachia', 'w')[:9]:>10}"
    for key in hosts:
        header += f"{'d|' + key[:7]:>11}"
    for key in hosts:
        header += f"{'wd|' + key[:7]:>11}"
    print(header)
    for design in results["designs"]:
        for sized in design["sizes"]:
            line = f"{design['id'][:43]:<44}{sized['requested_size']:>4}"
            got = sized["delivered_size"]
            line += f"{got if got is not None else 'n/a':>4}"
            evaluation = sized.get("evaluation")
            if evaluation is None:
                print(line + "  " + str(sized.get("evaluation_unavailable")))
                continue
            for key in targets:
                line += f"{_value(evaluation['per_target'][key]['coverage'], '{:.3f}'):>10}"
            reference = results["reference"]
            for field in ("selectivity_density", "weighted_selectivity_density"):
                for key in hosts:
                    pair = next(
                        p
                        for p in evaluation["target_host_pairs"]
                        if p["target"] == reference and p["host"] == key
                    )
                    line += f"{_value(pair[field], '{:.1f}'):>11}"
            print(line)
    print()
    print(f"Density columns are the reference strain ({results['reference']}) against each host:")
    print("d| from exact site counts (the Phase 2 pair figure), wd| from the weighted loads.")
    print()
    print(
        f"{'step':<60}{'cpus':>5}{'status':>9}{'seconds':>9}{'tree MiB':>10}"
        f"{'one MiB':>9}{'free':>10}"
    )
    rows = [("eval_tables count-kmers", results["eval_tables"])]
    for design in results["designs"]:
        for name, record in design.get("steps", {}).items():
            rows.append((f"{design['id']} {name}", record))
        for sized in design["sizes"]:
            rows.append((f"{design['id']} optimize n={sized['requested_size']}", sized["step"]))
    for name, record in rows:
        if not record or "wall_seconds" not in record:
            continue
        peak = record.get("peak_tree_rss_bytes")
        peak_text = "n/a" if peak is None else f"{peak / 2**20:.0f}"
        single = record.get("peak_single_process_rss_bytes")
        single_text = "n/a" if single is None else f"{single / 2**20:.0f}"
        free = f"{record.get('free_percent_before')}>{record.get('free_percent_after')}"
        cpus = record.get("cpus")
        print(
            f"{name[:59]:<60}{cpus if cpus is not None else 'n/a':>5}{record['status']:>9}"
            f"{record['wall_seconds']:>9.1f}{peak_text:>10}{single_text:>9}{free:>10}"
        )


# ----------------------------------------------------------------------
# The tables of the validation record
# ----------------------------------------------------------------------


def _pair_of(evaluation, target, host, field):
    """One pair figure of one evaluation, or None when the pair is not there."""
    if not evaluation:
        return None
    for pair in evaluation["target_host_pairs"]:
        if pair["target"] == target and pair["host"] == host:
            return pair.get(field)
    return None


def _measurement_of(evaluation, block, key, field):
    if not evaluation:
        return None
    return (evaluation.get(block) or {}).get(key, {}).get(field)


def _num(measurement, fmt="{:.2f}"):
    if measurement is None or measurement.get("value") is None:
        return "n/a"
    return fmt.format(measurement["value"])


def _spread_text(entry, fmt="{:.2f}"):
    """A C2 summary as `median [min, max]`, or why there is no figure.

    A seed whose panel binds the host nowhere carries the package's ceiling
    value rather than a measured ratio, and a median over a set holding one is
    not a median of measured ratios. Such a summary is reported as the count of
    those seeds instead of a number: the spread is not available, and printing
    one formed partly from the ceiling would read as though it were.
    """
    if entry is None or entry.get("unavailable"):
        return "n/a"
    zero = entry.get("seeds_with_zero_host_sites") or 0
    if zero:
        return f"ceiling in {zero}/{entry['n']}"
    return (
        f"{fmt.format(entry['median'])} [{fmt.format(entry['min'])}, "
        f"{fmt.format(entry['max'])}]"
    )


def _c2_blocks(sized):
    """The C2 summaries of one size, whichever pass measured each reference.

    The design's own hosts were summarised by the main run and the counts-only
    references by the counts-only pass. They are different blocks of the same
    file and are read as one mapping here, with nothing merged or recomputed.
    """
    out = {}
    for name in ("c2", "c2_counts_only"):
        block = sized.get(name)
        if block and block.get("summary"):
            out[name] = block["summary"]
    return out


def _c2_pair(sized, target, host, field):
    key = f"{target}|{host}"
    for summary in _c2_blocks(sized).values():
        entry = (summary.get(field) or {}).get(key)
        if entry is not None:
            return entry
    return None


def _c2_coverage(sized, target):
    for summary in _c2_blocks(sized).values():
        entry = (summary.get("coverage_by_target") or {}).get(target)
        if entry is not None:
            return entry
    return None


def _evaluations(sized):
    """The set-0 evaluations of one size, in the order they should be searched."""
    return [sized.get("evaluation"), sized.get("counts_only")]


def _find_pair(sized, target, host, field):
    for evaluation in _evaluations(sized):
        value = _pair_of(evaluation, target, host, field)
        if value is not None:
            return value
    return None


def _find_reference(sized, block, key, field):
    for evaluation in _evaluations(sized):
        value = _measurement_of(evaluation, block, key, field)
        if value is not None:
            return value
    return None


def print_tables(results):
    """Every table of the validation record, read from the results file.

    Nothing is recomputed from a genome here: each figure is the `Measurement`
    the Phase 2 function returned, or the spread an earlier pass formed over
    the C2 seeds. `n/a` is a figure that was not measured, and the JSON carries
    the reason beside it.
    """
    targets = results["targets"]
    reference = results["reference"]
    hosts = list(results["hosts"]) + list(results.get("counts_only_hosts") or [])
    genomes = results["genomes"]

    print("Table 1. The instance.")
    print(f"  reference target      {reference}")
    print(f"  targets               {', '.join(targets)}")
    print(f"  design hosts          {', '.join(results['hosts'])}")
    print(f"  counts-only hosts     {', '.join(results.get('counts_only_hosts') or []) or '-'}")
    print(f"  k                     {results['k']}")
    print(f"  sizes requested       {results['sizes_requested']}")
    print(f"  optimize seed         {results['optimize_seed']}")
    print(f"  C2 seeds              {results['c2_seeds']}")
    print(
        f"  polymerase, temp      {results['design_parameters']['polymerase']}, "
        f"{results['design_parameters']['reaction_temp']} C"
    )
    print(f"  git sha               {results['git_sha']}")
    print()
    print("Table 2. The genomes, as the manifest records them.")
    print(f"  {'key':<24}{'role':>8}{'bp':>14}{'records':>9}")
    for key in targets + hosts:
        genome = genomes[key]
        print(
            f"  {key:<24}{genome['role']:>8}{genome['length']:>14,}"
            f"{genome['records'] if genome['records'] is not None else 'n/a':>9}"
        )
    print()

    print("Table 3. Per-strain coverage of set 0, with the C2 median and range.")
    print("  Coverage is the union of per-primer reach windows. A reference read")
    print("  from an index and one scanned are not comparable (see the record).")
    header = f"  {'design':<42}{'req':>4}{'got':>4}"
    for key in targets:
        header += f"{key.replace('wolbachia', 'w')[:10]:>22}"
    print(header)
    for entry in results["designs"]:
        for sized in entry["sizes"]:
            line = f"  {entry['id'][:41]:<42}{sized['requested_size']:>4}"
            line += (
                f"{sized['delivered_size'] if sized['delivered_size'] is not None else 'n/a':>4}"
            )
            for key in targets:
                line += (
                    f"{_num(_find_reference(sized, 'per_target', key, 'coverage'), '{:.4f}'):>22}"
                )
            print(line)
            c2 = f"  {'  C2 median [min, max]':<42}{'':>8}"
            for key in targets:
                c2 += f"{_spread_text(_c2_coverage(sized, key), '{:.3f}'):>22}"
            print(c2)
    print()

    print("Table 4. The reference strain against each host: two different densities.")
    print("  exact  = the Phase 2 pair figure, from exact site counts per bp")
    print("  weight = derived from the occupancy-weighted loads, the quantity the")
    print("           optimizer reports as selectivity_density in occupancy mode")
    print("  zero   = C2 seeds with no host sites at all. The exact figure is then")
    print("           the package ceiling, not a measured ratio.")
    for host in hosts:
        print()
        print(f"  host {host} ({genomes[host]['length']:,} bp)")
        print(
            f"    {'design':<42}{'req':>4}{'exact':>10}{'C2 exact median [min, max]':>30}"
            f"{'zero':>6}{'weight':>10}{'C2 weight median [min, max]':>31}"
        )
        for entry in results["designs"]:
            for sized in entry["sizes"]:
                exact = _find_pair(sized, reference, host, "selectivity_density")
                weight = _find_pair(sized, reference, host, "weighted_selectivity_density")
                c2_exact = _c2_pair(sized, reference, host, "selectivity_density_by_pair")
                c2_weight = _c2_pair(sized, reference, host, "weighted_selectivity_density_by_pair")
                zero = "n/a" if c2_exact is None else c2_exact.get("seeds_with_zero_host_sites")
                print(
                    f"    {entry['id'][:41]:<42}{sized['requested_size']:>4}"
                    f"{_num(exact, '{:.3f}'):>10}{_spread_text(c2_exact, '{:.3f}'):>30}"
                    f"{zero:>6}{_num(weight, '{:.3f}'):>10}"
                    f"{_spread_text(c2_weight, '{:.3f}'):>31}"
                )
    print()

    print("Table 5. Host sites per Mb of set 0, exact counts.")
    print(f"  {'design':<42}{'req':>4}" + "".join(f"{h[:13]:>14}" for h in hosts))
    for entry in results["designs"]:
        for sized in entry["sizes"]:
            line = f"  {entry['id'][:41]:<42}{sized['requested_size']:>4}"
            for host in hosts:
                line += f"{_num(_find_reference(sized, 'per_host', host, 'sites_per_mb'), '{:.3f}'):>14}"
            print(line)
    print()

    print("Table 6. Delivered sets that are identical, by primer list.")
    seen = {}
    for entry in results["designs"]:
        for sized in entry["sizes"]:
            primers = sized.get("delivered_primers")
            if not primers:
                continue
            seen.setdefault(tuple(sorted(primers)), []).append(
                f"{entry['id']} n={sized['requested_size']}"
            )
    shared = {key: names for key, names in seen.items() if len(names) > 1}
    if not shared:
        print("  no two delivered sets are identical")
    for key, names in sorted(shared.items(), key=lambda item: -len(item[1])):
        print(f"  {len(names)} designs deliver the same {len(key)} primers:")
        for name in names:
            print(f"    {name}")
        print(f"    {', '.join(key)}")
    print()

    print("Table 7. Against the project's recorded figure for this pool.")
    recorded = results["recorded_reference"]
    print(
        f"  recorded: selectivity density {recorded['selectivity_density']} at coverage "
        f"{recorded['coverage']}, {recorded['panel_size']} primers"
    )
    print(f"  path:     {recorded['path']}")
    print(f"  pool:     {recorded['pool']}")
    print(f"  source:   {recorded['source']}")
    for entry in results["designs"]:
        comparison = entry.get("example_pool_comparison")
        if comparison:
            print(f"  {entry['id']} candidate pool against that pool: {comparison}")
        if entry["fg"] != [reference] or entry["bg"] != ["drosophila"]:
            continue
        for sized in entry["sizes"]:
            if sized["requested_size"] != recorded["panel_size"]:
                continue
            print(
                f"  this run, {entry['id']} n={sized['requested_size']}: "
                f"exact density {_num(_find_pair(sized, reference, 'drosophila', 'selectivity_density'), '{:.3f}')}, "
                f"weighted density {_num(_find_pair(sized, reference, 'drosophila', 'weighted_selectivity_density'), '{:.3f}')}, "
                f"coverage {_num(_find_reference(sized, 'per_target', reference, 'coverage'), '{:.4f}')}"
            )
    print()

    print("Table 8. What each step cost.")
    print_report(results)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) == 3 and argv[0] == "--worker":
        return worker_main(argv[1], argv[2])

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", default=os.path.join(REPO_ROOT, DEFAULT_MANIFEST))
    parser.add_argument(
        "--workdir",
        default=os.path.join(REPO_ROOT, DEFAULT_WORKDIR),
        help="Gitignored directory for tables, indexes and results.",
    )
    parser.add_argument("--output", help="Results JSON (default: WORKDIR/results.json).")
    parser.add_argument(
        "--hosts",
        default=",".join(DEFAULT_HOSTS),
        help="Host keys from the manifest, comma separated. Hosts not named are not opened.",
    )
    parser.add_argument(
        "--max-host-bp",
        type=int,
        default=1_000_000_000,
        help="Refuse a named host longer than this (default 1 Gb).",
    )
    parser.add_argument(
        "--counts-only-hosts",
        help="Host keys measured by k-mer counts alone, comma separated. They are "
        "evaluation references and nothing else: no design takes one as a "
        f"background and none is scanned for positions. A named host longer than "
        f"{COUNTS_ONLY_ABOVE_BP:,} bp is treated this way whichever flag named it.",
    )
    parser.add_argument(
        "--counts-only-pass",
        action="store_true",
        help="Add the counts-only references to an existing results file and exit. "
        "Builds no design and rebuilds no pool; the designs and step records "
        "already recorded are carried through.",
    )
    parser.add_argument(
        "--counts-only-tables-only",
        action="store_true",
        help="With --counts-only-pass, build the counts-only k-mer tables, record "
        "what they cost, and stop before evaluating anything. The count of a 3.3 Gb "
        "reference is the largest step of this stage and its peak is worth reading "
        "before more is started.",
    )
    parser.add_argument(
        "--counts-only-cpus",
        type=int,
        default=DESIGN_PARAMS["cpus"],
        help="Workers for the counts-only k-mer count "
        f"(default {DESIGN_PARAMS['cpus']}, as every other step).",
    )
    parser.add_argument(
        "--counts-route-check",
        metavar="DESIGN_ID",
        help="Measure one host of this design's set both ways and compare the site "
        "counts: the counts route against the figure already recorded through "
        "positions. Appends the comparison to the results file and exits non-zero "
        "if the two disagree.",
    )
    parser.add_argument(
        "--counts-route-check-host",
        default="drosophila",
        help="Host for --counts-route-check (default drosophila).",
    )
    parser.add_argument(
        "--counts-route-check-size",
        type=int,
        default=12,
        help="Requested panel size for --counts-route-check (default 12).",
    )
    parser.add_argument(
        "--tables",
        action="store_true",
        help="Print the tables of the validation record from an existing results "
        "file and exit. Reads; measures nothing.",
    )
    parser.add_argument("--c1-host", help="Host for the leave-one-out designs (default: largest).")
    parser.add_argument("--sizes", default=",".join(str(size) for size in DEFAULT_SIZES))
    parser.add_argument("--only", help=f"Groups to run, comma separated, from {','.join(GROUPS)}.")
    parser.add_argument("--match", help="Only designs whose id contains this text.")
    parser.add_argument("--seeds", type=int, default=20, help="C2 panels per design and size.")
    parser.add_argument(
        "--c2-scan-offdesign-hosts",
        action="store_true",
        help="Scan hosts outside the design for C2 panels too. Without it their sites "
        "come from k-mer counts and their coverage is unavailable.",
    )
    parser.add_argument(
        "--pooled-cpus",
        type=int,
        help="Workers for steps 1 to 3 of a design with more than one foreground genome "
        f"(default {DESIGN_PARAMS['cpus']}, as every other design). optimize is not affected.",
    )
    parser.add_argument(
        "--cpus-check",
        metavar="DESIGN_ID",
        help="Rebuild this design's pool at --pooled-cpus in a separate directory, compare "
        "step2_df.csv and step3_df.csv byte for byte with the design's own, and exit.",
    )
    parser.add_argument("--retry", action="store_true", help="Re-run stopped or failed steps.")
    parser.add_argument("--rss-limit-gb", type=float, default=9.0)
    parser.add_argument(
        "--halt-peak-gb",
        type=float,
        default=8.0,
        help="Halt the whole run after a step whose peak passed this figure.",
    )
    parser.add_argument("--stop-free-percent", type=int, default=20)
    parser.add_argument("--start-free-percent", type=int, default=60)
    parser.add_argument("--start-wait-seconds", type=int, default=600)
    parser.add_argument("--min-free-disk-gb", type=float, default=15.0)
    parser.add_argument("--sample-seconds", type=float, default=1.0)
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
