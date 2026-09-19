#!/bin/bash
# One panel size, two frontiers, through the real `optimize` command.
#
# The only thing that differs between the two runs is step3_df.csv: the
# shortlist `optimize` reads today, against the whole eligible inventory a
# frontier refill would reach. Same seed, same params, same position files.
set -euo pipefail

SP="$1"
SIZE="$2"
SEED="${3:-20260919}"
WORK="$SP/frontier/work"
OUTDIR="$SP/frontier/results/n${SIZE}_s${SEED}"
mkdir -p "$OUTDIR"

cp "$WORK/step3_df.csv" "$SP/frontier/step3_narrow.csv"

for MODE in narrow wide; do
  cp "$SP/frontier/step3_$MODE.csv" "$WORK/step3_df.csv"
  rm -f "$WORK/step4_improved_df.csv" "$WORK/step4_improved_df_summary.json"
  echo "### $MODE  n=$SIZE  pool=$(( $(wc -l < "$WORK/step3_df.csv") - 1 ))"
  START=$(python -c 'import time; print(time.monotonic())')
  neoswga optimize -j "$SP/frontier/params.json" -n "$SIZE" --seed "$SEED" \
    > "$OUTDIR/$MODE.log" 2>&1 || { echo "FAILED: see $OUTDIR/$MODE.log"; tail -20 "$OUTDIR/$MODE.log"; exit 1; }
  END=$(python -c 'import time; print(time.monotonic())')
  python -c "print(f'    {$END - $START:.1f}s')"
  cp "$WORK/step4_improved_df.csv" "$OUTDIR/$MODE.step4.csv"
  cp "$WORK/step4_improved_df_summary.json" "$OUTDIR/$MODE.summary.json"
  python -c "print(f'    seconds: {$END - $START:.1f}')" > "$OUTDIR/$MODE.seconds"
done

cp "$SP/frontier/step3_narrow.csv" "$WORK/step3_df.csv"
echo "done: $OUTDIR"
