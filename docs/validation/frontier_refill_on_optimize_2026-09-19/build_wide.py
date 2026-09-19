"""Write a step3_df.csv holding the WHOLE eligible inventory, not the shortlist.

`optimize` opens its frontier at the size of the list it reads, so swapping
this file is the same lever a refill would pull, exercised through the real
command rather than through a shim around the optimizer.

Only the `primer` column reaches the optimizer; the rest are carried so the
file is the shape every other reader expects. The ORDER is the inventory's own
`search_rank`, which is what a widened frontier would hand the search -- not
the step-2 ranking, which stops at `max_primer`.

The condition and QC policy are read out of the inventory rather than
re-derived from params.json. Re-deriving them means rebuilding the reaction
fingerprint, and a fingerprint that misses returns an empty set silently, which
is exactly the failure this measurement must not make.
"""

import os
import sqlite3
import sys

import pandas as pd

from neoswga.core.candidate_inventory import CandidateInventory

DATA_DIR = sys.argv[1]
OUT = sys.argv[2]

db_path = os.path.join(DATA_DIR, "candidate_inventory.sqlite")

# The hard-QC assessment, not the `stage3:carry_forward` one, which is the
# shortlist by another name and would reproduce the narrow frontier.
with sqlite3.connect(db_path) as conn:
    keys = conn.execute(
        "SELECT condition_id, policy_version, COUNT(*) FROM assessments "
        "WHERE passed = 1 AND policy_version LIKE 'qc-%' "
        "GROUP BY 1, 2 ORDER BY 3 DESC"
    ).fetchall()
if not keys:
    raise SystemExit("no hard-QC assessments in the inventory")
condition_id, policy_version, n_passed = keys[0]
print(f"condition: {condition_id}")
print(f"policy:    {policy_version}  ({n_passed:,} passed)")

narrow = pd.read_csv(os.path.join(DATA_DIR, "step3_df.csv"))
lengths = sorted({len(p) for p in narrow["primer"].astype(str)})

inv = CandidateInventory(db_path)
eligible = list(inv.iter_eligible(condition_id, lengths, policy_version=policy_version))
print(f"eligible: {len(eligible):,}   shortlist: {len(narrow):,}")
if len(eligible) <= len(narrow):
    raise SystemExit("the wide universe is no wider than the shortlist; nothing to measure")

missing = set(narrow["primer"].astype(str)) - set(eligible)
print(f"shortlist primers absent from the eligible set: {len(missing):,}")

step2 = pd.read_csv(os.path.join(DATA_DIR, "step2_df.csv"))
cols = [c for c in ("fg_count", "bg_count", "gini", "ratio") if c in step2.columns]
lookup = step2.set_index("primer")[cols]

wide = pd.DataFrame({"primer": eligible})
wide["step2_rank"] = range(len(wide))
wide = wide.join(lookup, on="primer")
wide = wide[["primer", "step2_rank", "ratio", "gini", "fg_count", "bg_count"]]
wide.to_csv(OUT, index=False)
print(f"wrote {OUT}: {len(wide):,} rows, {int(wide['fg_count'].notna().sum()):,} carry step-2 measurements")
