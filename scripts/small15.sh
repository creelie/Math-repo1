#!/usr/bin/env bash
# Run the K = 15 small-prime gates (17 <= p <= 31) and keep their summaries.
# usage: scripts/small15.sh WORKERS p1 p2 ...
#   WORK    scratch directory (default ./work); KEEP=1 keeps the row files
# Each gate is resumable: rerunning reuses the finished chunks in WORK/small_<p>.
set -u
[ $# -ge 2 ] || { echo "usage: $0 WORKERS p1 p2 ..." >&2; exit 2; }
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=${WORK:-$ROOT/work}
KEEP=${KEEP:-0}
W=$1; shift
[ -x "$ROOT/bin/cascade_k15" ] || { echo "missing bin/cascade_k15 (run make)" >&2; exit 2; }
mkdir -p "$WORK" "$ROOT/results"
for p in "$@"; do
  out=$WORK/small_$p
  rm -f "$out/SUMMARY.json" "$out/MANIFEST_SHA256.json"
  start=$(date +%s)
  python3 "$ROOT/src/small_gate.py" "$ROOT/bin/cascade_k15" "$p" "$W" "$out" >> "$WORK/s_$p.log" 2>&1
  rc=$?
  echo "p=$p rc=$rc wall=$(( $(date +%s) - start ))s $(tail -1 "$WORK/s_$p.log")" | tee -a "$WORK/SWEEP.txt"
  [ -f "$out/SUMMARY.json" ] || continue
  cp "$out/SUMMARY.json" "$ROOT/results/summary_$p.json"
  cp "$out/MANIFEST_SHA256.json" "$ROOT/results/manifest_$p.json"
  [ "$KEEP" = 1 ] || rm -f "$out"/rows_*.txt "$out"/filt_*.out "$out"/filt_*.json
done
