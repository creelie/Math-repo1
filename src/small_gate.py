#!/usr/bin/env python3
"""Prime gate for K = 15 at the primes 17 <= p <= 31, without the cover generator.

usage: small_gate.py CASCADE P WORKERS OUTDIR

For p < 32 a speed class c covers a time T at level one only when
16 * dist(T c, p) < p, that is dist(T c, p) = 1, so c covers exactly the
times T = +-c^{-1} (mod p). A level-one cover of the (p-1)/2 folded times
therefore contains every class, and the 15 - (p-1)/2 remaining speeds are
arbitrary classes. The script checks this lemma on the actual covering
relation, lists every such multiset, keeps one representative of each orbit
under multiplication by units modulo p (which maps the lift tree of a row
onto that of its image), and runs the exact binary lift cascade on the
representatives. The gate closes when no row keeps an improper lift at the
cascade's deepest level.

The representatives are split into chunks of at most CHUNK rows, independent
of WORKERS. A finished chunk leaves OUTDIR/filt_<i>.json; rerunning with the
same OUTDIR checks that the chunk files are unchanged and runs only the
chunks without it, so an interrupted gate resumes where it stopped.

Writes OUTDIR/SUMMARY.json and OUTDIR/MANIFEST_SHA256.json; exit 0 if closed,
2 otherwise.
"""
import hashlib, itertools, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

K = 15
CHUNK = 4
if len(sys.argv) != 5:
    print(__doc__); sys.exit(1)
CASCADE, P, WORKERS, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
if P < 17 or P > 31 or any(P % q == 0 for q in range(2, P)):
    sys.exit("P must be a prime with 17 <= P <= 31")
NC = (P - 1) // 2
FREE = K - NC
os.makedirs(OUT, exist_ok=True)
t_start = time.time()

def dist(x):
    x %= P
    return min(x, P - x)

# the covering relation at level one, and the lemma that makes every class necessary
for c in range(1, NC + 1):
    for T in range(1, NC + 1):
        covers = (K + 1) * dist(T * c) < P
        if covers != (dist(T * c) == 1):
            sys.exit(f"lemma fails at c={c} T={T}")
for T in range(1, NC + 1):
    if sum((K + 1) * dist(T * c) < P for c in range(1, NC + 1)) != 1:
        sys.exit(f"time {T} is not covered by exactly one class")

FOLD = [min(x, P - x) for x in range(P)]
def canon(row):
    return min(tuple(sorted(FOLD[g * c % P] for c in row)) for g in range(1, P))

orbits = {}
n_rows = 0
for extra in itertools.combinations_with_replacement(range(1, NC + 1), FREE):
    row = tuple(sorted(list(range(1, NC + 1)) + list(extra)))
    n_rows += 1
    key = canon(row)
    orbits[key] = orbits.get(key, 0) + 1
reps = sorted(orbits)

NCH = -(-len(reps) // CHUNK)
files = []
for i in range(NCH):
    fn = f"{OUT}/rows_{i}.txt"
    text = "".join("".join(f"{c:02x}" for c in r) + "\n" for r in reps[i::NCH])
    if os.path.exists(fn):
        if open(fn).read() != text: sys.exit(f"{fn} differs from the rows of this gate; use a fresh OUTDIR")
    else:
        open(fn + ".tmp", "w").write(text); os.replace(fn + ".tmp", fn)
    files.append(fn)

def parse(fn, out, err, secs):
    m = re.search(r"TOTAL rows=(\d+) survivors=(\d+) shiftkilled=(\d+) maxlevel=(\d+)", open(err).read())
    if not m: sys.exit(f"cascade stats without maxlevel on {fn}")
    return dict(file=os.path.basename(fn), rows=int(m[1]), survivors=int(m[2]), shift=int(m[3]),
                level=int(m[4]), secs=secs, out=out, stats=err)

def names(fn):
    base = fn.replace("rows_", "filt_")[:-len(".txt")]
    return base + ".out", base + ".stats", base + ".json"

def run(fn):
    out, err, done = names(fn)
    if os.path.exists(done):
        return parse(fn, out, err, json.load(open(done))["secs"]) | {"resumed": True}
    t0 = time.time()
    r = subprocess.run([CASCADE, str(P), "filter", fn], capture_output=True, text=True)
    if r.returncode != 0: sys.exit(f"cascade failed on {fn}: {r.stderr[-300:]}")
    for path, text in ((out, r.stdout), (err, r.stderr)):
        open(path + ".tmp", "w").write(text); os.replace(path + ".tmp", path)
    secs = round(time.time() - t0, 1)
    x = parse(fn, out, err, secs)
    open(done + ".tmp", "w").write(json.dumps(dict(secs=secs))); os.replace(done + ".tmp", done)
    print(f"chunk {os.path.basename(fn)} rows={x['rows']} survivors={x['survivors']} secs={secs}", flush=True)
    return x | {"resumed": False}

with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    res = list(ex.map(run, files))
resumed = sum(x.pop("resumed") for x in res)
if sum(x["rows"] for x in res) != len(reps): sys.exit("cascade row count mismatch")
levels = {x["level"] for x in res}
if len(levels) != 1: sys.exit("mixed lift levels")
LEVEL = levels.pop()
JMAX = LEVEL.bit_length() - 1

tally, persistent = {}, []
for x in res:
    for line in open(x["out"]):
        if not line.startswith("SURVIVOR"): continue
        if "PERSISTENT" in line:
            persistent.append(line.strip()); continue
        for lv, cnt in re.findall(r" l(\d+)=(\d+)", line)[1:]:
            if cnt == "0":
                tally[f"died_l{lv}"] = tally.get(f"died_l{lv}", 0) + 1; break
closed = not persistent
shift = sum(x["shift"] for x in res)
claim = (f"Assuming LRC(m) for all m < {K} (arXiv:2609.02604v2), every primitive counterexample v to LRC({K}) "
         f"({K} speeds, {K+1} runners) with pairwise distinct coordinates has {P} | v_1...v_{K}. "
         f"Every level-one cover modulo {P} contains all {NC} classes; each of its lifts at level {LEVEL} is "
         f"({K},{P},{LEVEL})-proper (witness on the 1/({LEVEL}*{P}) grid or the gcd clause) or satisfies the "
         f"power-of-two shift lemma (some d = 2^j, 2 <= j <= {JMAX}); {shift} lifts of the orbit representatives use "
         f"the lemma, so this is the divisibility of Lemma 2.2(ii) of arXiv:2609.02604v2.")
summary = dict(
    gate=dict(k=K, p=P, lift_level=LEVEL, method="small prime: every class is needed, unit-orbit representatives"),
    claim=claim if closed else None, status="GATE_CLOSED" if closed else "NOT_CLOSED",
    covers=dict(classes=NC, free_speeds=FREE, multisets=n_rows, orbit_representatives=len(reps)),
    cascade=dict(rows=len(reps), survivor_rows=sum(x["survivors"] for x in res), lifts_settled_by_shift_lemma=shift,
                 persistent_rows=len(persistent), **tally),
    persistent=persistent[:50],
    jobs=[{k: v for k, v in x.items() if k not in ("out", "stats")} for x in res],
    chunks=dict(count=NCH, rows_per_chunk_at_most=CHUNK, reused_from_earlier_runs=resumed),
    wall_secs_last_run=round(time.time() - t_start, 1),
    tools=dict(cascade=os.path.basename(CASCADE), script=os.path.basename(__file__)),
)
json.dump(summary, open(f"{OUT}/SUMMARY.json", "w"), indent=2)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
manifest = {os.path.basename(f): sha256(f) for f in files + [x["out"] for x in res] + [x["stats"] for x in res]}
manifest[os.path.basename(CASCADE)] = sha256(CASCADE)
manifest[os.path.basename(__file__)] = sha256(__file__)
json.dump(manifest, open(f"{OUT}/MANIFEST_SHA256.json", "w"), indent=2)
print(f"{summary['status']} p={P} level={LEVEL} orbits={len(reps)} multisets={n_rows} persistent={len(persistent)} "
      f"chunks={NCH} reused={resumed} wall={summary['wall_secs_last_run']}s")
sys.exit(0 if closed else 2)
