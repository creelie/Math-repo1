#!/usr/bin/env bash
# Run the K = 15 prime gates in order and keep their summaries.
# usage: scripts/sweep15.sh WORKERS p1 p2 ...
#   WORK    scratch directory for rows and cascade outputs (default ./work)
#   KEEP=1  keep the bulky rows and cascade outputs (default: delete them
#           after copying SUMMARY.json and MANIFEST_SHA256.json)
# Build first with `make`. Each gate is resumable: rerunning reuses finished
# jobs in WORK/out_<p>.
set -u
[ $# -ge 2 ] || { echo "usage: $0 WORKERS p1 p2 ..." >&2; exit 2; }
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=${WORK:-$ROOT/work}
KEEP=${KEEP:-0}
W=$1; shift
for tool in python3 "$ROOT/bin/cascade_k15" "$ROOT/bin/bgk15_nw2"; do
  command -v "$tool" >/dev/null 2>&1 || [ -x "$tool" ] || { echo "missing $tool (run make)" >&2; exit 2; }
done
mkdir -p "$WORK" "$ROOT/results"
for p in "$@"; do
  n=$(( (p - 1) / 2 )); nw=$(( (n + 63) / 64 )); [ $nw -lt 2 ] && nw=2
  eng=$ROOT/bin/bgk15_nw$nw
  [ -x "$eng" ] || { echo "no engine for p=$p (needs NW=$nw; add it to the Makefile)" >&2; exit 2; }
  out=$WORK/out_$p
  start=$(date +%s)
  python3 "$ROOT/src/gate_k15.py" "$eng" "$ROOT/bin/cascade_k15" /bin/true "$p" "$W" "$out" --no-kill \
    >> "$WORK/g_$p.log" 2>&1
  rc=$?
  verdict=$(grep -o -E "GATE_CLOSED|NOT_CLOSED|FAIL.*" "$WORK/g_$p.log" | tail -1)
  echo "p=$p rc=$rc wall=$(( $(date +%s) - start ))s verdict=$verdict" | tee -a "$WORK/SWEEP.txt"
  cp "$out/SUMMARY.json" "$ROOT/results/summary_$p.json" 2>/dev/null
  cp "$out/MANIFEST_SHA256.json" "$ROOT/results/manifest_$p.json" 2>/dev/null
  [ "$KEEP" = 1 ] || rm -rf "$out/rows_ir" "$out/rows_km1" "$out/filt"
done
