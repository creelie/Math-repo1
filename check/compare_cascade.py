#!/usr/bin/env python3
"""Compare the cascade with the independent checker row by row.

usage: compare_cascade.py P ROWFILE [CASCADE] [CHECKER]

ROWFILE holds rows in the cascade's input format (two or four hex digits per
class). Both programs walk the lift tree of every row to the cascade's depth
(its `maxlevel`); the script requires identical counts of unsettled lifts at
levels 4, 8, ..., identical shift-lemma counts and the same verdict. The
level-2 count may differ: the cascade counts level-2 lifts before the gcd
clause, the checker after it.
"""
import re, subprocess, sys

if len(sys.argv) < 3:
    print(__doc__); sys.exit(2)
P, rowfile = int(sys.argv[1]), sys.argv[2]
cascade = sys.argv[3] if len(sys.argv) > 3 else "bin/cascade_k15"
checker = sys.argv[4] if len(sys.argv) > 4 else "bin/lift_check_sym"
K = 15

r = subprocess.run([cascade, str(P), "filter", rowfile], capture_output=True, text=True)
if r.returncode != 0: sys.exit(f"cascade failed: {r.stderr[-300:]}")
m = re.search(r"maxlevel=(\d+)", r.stderr)
if not m: sys.exit("cascade without maxlevel")
level = int(m[1]); depth = level.bit_length() - 1
casc = {}
for line in r.stdout.splitlines():
    if not line.startswith("SURVIVOR"): continue
    base = tuple(int(x) for x in line.split()[2:2 + K])
    casc[base] = line

def counts(line):
    c = {int(a): int(b) for a, b in re.findall(r" l(\d+)=(\d+)", line)}
    s = int(re.search(r"shift=(\d+)", line)[1])
    return c, s, "PERSISTENT" in line

ok = bad = 0
for raw in open(rowfile):
    raw = raw.strip()
    if not raw: continue
    w = 4 if len(raw) == 4 * K else 2
    base = tuple(int(raw[i:i + w], 16) for i in range(0, len(raw), w))
    out = subprocess.run([checker, str(depth), str(P)] + [str(x) for x in base],
                         capture_output=True, text=True).stdout
    c2, s2, p2 = counts(out)
    if base not in casc:          # the cascade prints only rows with level-2 improper lifts
        good = all(v == 0 for k, v in c2.items() if k >= 4) and s2 == 0 and not p2
    else:
        c1, s1, p1 = counts(casc[base])
        lv = [1 << j for j in range(2, depth + 1)]
        good = ([c1.get(k, 0) for k in lv] == [c2.get(k, 0) for k in lv] and s1 == s2 and p1 == p2)
    if good: ok += 1
    else:
        bad += 1
        print("MISMATCH", base, casc.get(base, "(no cascade line)"), "|", out.strip())
print(f"compare p={P} depth={level}: {ok} rows agree, {bad} differ")
sys.exit(1 if bad else 0)
