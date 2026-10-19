#!/usr/bin/env python3
"""Complete K=15 prime-gate pipeline for one prime p (LRC(15), sixteen runners).
Ported from Allikvere's gate_k14.py (10.5281/zenodo.22667683, CC BY 4.0).
Stages (all exact; any failure or surviving improper lift halts with exit 2):
  1. precondition  ENGINE p km1low 13  -> canonical covers on <= 13 classes.
                   Quota-pruned, class-1 rooted (some class of a cover has
                   enough private times, and a unit moves it to class 1).
                   No such cover (tau(p) >= 14) licenses the decomposition of
                   the reducible branch into exactly-14-class (hence minimal)
                   covers; otherwise the reducible branch is the single
                   general km1all job.
  2. irredundant   ENGINE p irsubrootrawc r s rows  for every (root, subroot):
                   orbit-canonical irredundant 15-class covers (branch 0).
  3. reducible     ENGINE p km1root r rows  for every km1 root (or km1all):
                   canonical 14-class covers; the cascade extends each by all
                   NC classes (repeated-class and redundant-support rows).
  4. cascade       CASCADE p filter rows (irredundant) / filterext rows (14-covers):
                   exact binary lift cascade l=2..128 with per-row exact fibers;
                   a lift is settled by a witness, the gcd clause (>= K-1 even
                   coordinates) or the power-of-two shift lemma (cascade
                   source, shift_proper; third alternative of Lemma 2.2(ii) of
                   arXiv:2609.02604v2).
  5. orbits        PERSISTENT survivor rows (a lift left at the deepest level,
                   128 by default; printed as maxlevel by the cascade)
                   deduplicated by unit orbit.  The gate closes only if there
                   are none; K+1 = 16 is a power of two, so there is no odd
                   terminal level and --no-kill is the K=15 setting.
  6. kill          (K=14 only) level-15 fiber split 3x5.
  7. summary       SUMMARY.json + MANIFEST_SHA256.json in OUTDIR.
usage: gate_k15.py ENGINE CASCADE KILL P WORKERS OUTDIR [--shard i/N] [--no-kill]
Sharding splits stages 2-4 (jobs with index % N == i); a sharded run writes
SHARD_SUMMARY.json and skips stages 5-7.
Resumable: finished job outputs are reused.
"""
import glob, hashlib, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

K = 15
args = [a for a in sys.argv[1:] if not a.startswith("--")]
flags = [a for a in sys.argv[1:] if a.startswith("--")]
if len(args) < 6:
    print(__doc__); sys.exit(1)
ENGINE, CASCADE, KILL, P, WORKERS, OUT = args[0], args[1], args[2], int(args[3]), int(args[4]), args[5]
shard_i, shard_n = 0, 1
for f in flags:
    if f.startswith("--shard"):
        shard_i, shard_n = map(int, f.split("=")[1].split("/"))
NO_KILL = "--no-kill" in flags
NC = (P - 1) // 2
os.makedirs(OUT, exist_ok=True)
for d in ("rows_ir", "rows_km1", "filt", "kills"):
    os.makedirs(f"{OUT}/{d}", exist_ok=True)
LOG = open(f"{OUT}/gate_{P}_shard{shard_i}of{shard_n}.log", "a")

def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def run(cmd, timeout=None):
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        return out, f"TIMEOUT after {timeout}s", 124, time.time() - t0
    return r.stdout, r.stderr, r.returncode, time.time() - t0

KILL_TIMEOUT = float(os.environ.get("KILL_TIMEOUT", "1800")) or None   # seconds per orbit kill (0 = no limit); campaign kills take < 1 s

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def fail(msg):
    log("FAIL " + msg); sys.exit(2)

def lift_level(stats, name):
    # deepest binary lift level of the cascade binary that produced `stats`
    m = re.search(r"maxlevel=(\d+)", stats)
    if not m: fail(f"{name}: cascade stats without maxlevel (old cascade binary?)")
    return int(m[1])

timing = {}
log(f"START p={P} K={K} shard={shard_i}/{shard_n} engine={ENGINE} cascade={CASCADE} kill={KILL}")
info, _, rc, _ = run([ENGINE, str(P), "info"])
log(info.strip())
if rc != 0 or f"K={K} " not in info: fail("engine info")

# ---------- 1. precondition (shard 0 only; others trust the merge step) ----------
pre = {}
if shard_i == 0:
    pf = f"{OUT}/km1low_{P}.out"
    if os.path.exists(pf) and "canonical=" in open(pf).read():
        out = open(pf).read()
    else:
        out, err, rc, secs = run([ENGINE, str(P), "km1low", str(K - 2), f"{OUT}/km1low_{P}.txt"])
        if rc != 0: fail("km1low " + err)
        open(pf, "w").write(out)
    m = re.search(r"km1low smax=(\d+) canonical=(\d+) nodes=(\d+) .*secs=([\d.]+)", out)
    if not m or int(m[1]) != K - 2: fail("km1low parse")
    pre = dict(smax=int(m[1]), covers_upto_smax=int(m[2]), nodes=int(m[3]), secs=float(m[4]))
    timing["precondition_secs"] = pre["secs"]
    log(f"PRECONDITION tau(p)>={K-1}: covers on <= {K-2} classes = {pre['covers_upto_smax']} ({pre['secs']:.1f}s)")
    # tau(p) >= K-1 licenses the decomposition: every (K-1)-class cover is then
    # minimal, so the irredundant km1root enumeration is the complete reducible
    # branch. If tau(p) <= K-2 (small p) that fails - non-irredundant (K-1)-class
    # covers exist and neither km1root (irredundant only) nor km1low (<= K-2
    # classes) produces them - so the general variant runs km1all (paper '12all'
    # for K=13): canonical (K-1)-multiset covers on quota-passing supports of
    # <= K-1 classes containing class 1, class-1 multiplicity one, quota q_{K-1}.
    # Complete for the one-class extensions by the quota lemma with the
    # multiplicity refinement (some private-rich class has multiplicity <= 2;
    # see README). Single serial job, affordable only at small p.
    pre["variant"] = "decomposition" if pre["covers_upto_smax"] == 0 else "general (km1all)"
    if pre["covers_upto_smax"] != 0:
        log(f"GENERAL REDUCIBLE VARIANT: tau(p) <= {K-2}, km1root replaced by the complete km1all family")
GENERAL = pre.get("variant", "").startswith("general")

# ---------- 2+3. job lists ----------
def lines(a):
    return run([ENGINE, str(P)] + a)[0].splitlines()
roots = [int(l.split()[0]) for l in lines(["irroots"])[1:]]
ir_jobs = []
for r in roots:
    subs = [int(l.split()[0]) for l in lines(["irsubroots", str(r)])[1:]]
    ir_jobs += [(r, s) for s in subs]
km1_roots = [int(l.split()[0]) for l in lines(["km1roots"])[1:]]
ir_jobs_all, km1_all = len(ir_jobs), len(km1_roots)
ir_jobs = [j for i, j in enumerate(ir_jobs) if i % shard_n == shard_i]
km1_jobs = [r for i, r in enumerate(km1_roots) if i % shard_n == shard_i]
if GENERAL:
    if shard_n > 1: fail("general reducible variant is not sharded: run this small gate unsharded")
    km1_jobs = ["all"]
log(f"jobs: irredundant {len(ir_jobs)}/{ir_jobs_all} (roots={len(roots)}), " +
    ("km1all 1 (general variant)" if GENERAL else f"km1root {len(km1_jobs)}/{km1_all}"))

# ---------- 2-4. generate + cascade per job ----------
def cascade(mode, rows, name):
    outf, statf = f"{OUT}/filt/{name}.out", f"{OUT}/filt/{name}.stats"
    if os.path.exists(statf) and "TOTAL rows=" in open(statf).read():
        return open(statf).read(), 0.0
    out, err, rc, secs = run([CASCADE, str(P), mode, rows])
    if rc != 0: fail(f"cascade {name} rc={rc} {err[-500:]}")
    open(outf, "w").write(out); open(statf, "w").write(err)
    return err, secs

def ir_job(job):
    r, s = job
    name = f"ir_{r}_{s}"; rows = f"{OUT}/rows_ir/{name}.txt"; lf = rows + ".log"
    if os.path.exists(lf) and "raw_irred" in open(lf).read():
        gout, gsecs = open(lf).read(), 0.0
    else:
        gout, err, rc, gsecs = run([ENGINE, str(P), "irsubrootrawc", str(r), str(s), rows])
        if rc != 0 or "raw_irred" not in gout: fail(f"gen {name} {err[-300:]}")
        open(lf, "w").write(gout)
    m = re.search(rf"raw_irred{K}=(\d+) rawscan=(\d+) nodes=(\d+) .*secs=([\d.]+)", gout)
    stats, csecs = cascade("filter", rows, name)
    m2 = re.search(r"TOTAL rows=(\d+) survivors=(\d+) shiftkilled=(\d+)", stats)
    if not m2: fail(f"{name}: cascade stats without shiftkilled (old cascade binary?)")
    if int(m2[1]) != int(m[1]): fail(f"{name}: cascade rows {m2[1]} != generated {m[1]}")
    return dict(name=name, rows=int(m[1]), leaves=int(m[2]), nodes=int(m[3]), gen_secs=float(m[4]),
                gen_wall=gsecs, survivors=int(m2[2]), shift=int(m2[3]), casc_wall=csecs, level=lift_level(stats, name))

def km1_job(r):
    # canonical 13-class covers rooted at class 1 whose second class is root r;
    # one orbit can appear under several roots (up to 13 rooted images), so the
    # union is deduplicated before the extension cascade.
    name = f"km1_{r}"; rows = f"{OUT}/rows_km1/{name}.txt"; lf = rows + ".log"
    cmd = [ENGINE, str(P), "km1all", rows] if r == "all" else [ENGINE, str(P), "km1root", str(r), rows]
    if os.path.exists(lf) and "canonical" in open(lf).read():
        gout, gsecs = open(lf).read(), 0.0
    else:
        gout, err, rc, gsecs = run(cmd)
        if rc != 0 or "canonical" not in gout: fail(f"gen {name} {err[-300:]}")
        open(lf, "w").write(gout)
    m = re.search(r"canonical\d+=(\d+) nodes=(\d+) leaves=(\d+) secs=([\d.]+)", gout)
    if r == "all" and int(m[1]) < pre["covers_upto_smax"]: fail(f"km1all {m[1]} covers < km1low {pre['covers_upto_smax']}")
    return dict(name=name, covers=int(m[1]), nodes=int(m[2]), gen_secs=float(m[4]), gen_wall=gsecs)

def ext_job(chunk):
    # chunk = (index, rows file with deduplicated 13-class covers)
    idx, rows = chunk
    n = sum(1 for _ in open(rows))
    stats, csecs = cascade("filterext", rows, f"ext_{idx}")
    m2 = re.search(r"TOTAL rows=(\d+) survivors=(\d+) shiftkilled=(\d+)", stats)
    if not m2: fail(f"ext_{idx}: cascade stats without shiftkilled (old cascade binary?)")
    if int(m2[1]) != n * NC: fail(f"ext_{idx}: extension rows {m2[1]} != {n}*{NC}")
    return dict(name=f"ext_{idx}", covers=n, ext_rows=int(m2[1]), survivors=int(m2[2]), shift=int(m2[3]), casc_wall=csecs,
                level=lift_level(stats, f"ext_{idx}"))

t0 = time.time()
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    ir_res = list(ex.map(ir_job, ir_jobs))
    km1_res = list(ex.map(km1_job, km1_jobs))
    covers = set()
    for x in km1_res:
        covers.update(l.strip() for l in open(f"{OUT}/rows_km1/{x['name']}.txt") if l.strip())
    if GENERAL:      # consistency: the km1low rows (<= K-2 classes) must all be in the km1all family
        low = set(l.strip() for l in open(f"{OUT}/km1low_{P}.txt") if l.strip())
        if len(low) != pre["covers_upto_smax"] or not low <= covers:
            fail(f"km1low rows {len(low)} (reported {pre['covers_upto_smax']}) not contained in km1all family of {len(covers)}")
    covers = sorted(covers)
    chunks = []
    nchunk = max(1, min(WORKERS * 2, len(covers)))
    for idx in range(nchunk):
        part = covers[idx::nchunk]
        if not part: continue
        fn = f"{OUT}/rows_km1/dedup_{idx}.txt"
        open(fn, "w").write("\n".join(part) + "\n")
        chunks.append((idx, fn))
    ext_res = list(ex.map(ext_job, chunks))
timing["gen_cascade_wall_secs"] = time.time() - t0
levels = {x["level"] for x in ir_res + ext_res}
if len(levels) != 1: fail(f"cascade outputs from binaries with different lift levels {sorted(levels)} (resumed with a new binary?)")
LEVEL = levels.pop()
JMAX = LEVEL.bit_length() - 1
timing["ir_gen_cpu_secs"] = sum(x["gen_secs"] for x in ir_res)
timing["ir_cascade_wall_secs"] = sum(x["casc_wall"] for x in ir_res)
timing["km1_gen_cpu_secs"] = sum(x["gen_secs"] for x in km1_res)
timing["km1_cascade_wall_secs"] = sum(x["casc_wall"] for x in ext_res)
base = dict(irredundant_jobs=len(ir_res), irredundant_rows=sum(x["rows"] for x in ir_res),
            irredundant_nodes=sum(x["nodes"] for x in ir_res),
            km1_jobs=len(km1_res), km1_rooted_covers=sum(x["covers"] for x in km1_res),
            km1_covers_dedup=len(covers), reducible_variant=pre.get("variant", "decomposition (shard>0)"),
            km1_extension_rows=sum(x["ext_rows"] for x in ext_res), extensions_per_cover=NC,
            level2_survivor_rows_ir=sum(x["survivors"] for x in ir_res),
            level2_survivor_rows_ext=sum(x["survivors"] for x in ext_res),
            lifts_settled_by_shift_lemma=sum(x["shift"] for x in ir_res) + sum(x["shift"] for x in ext_res))
if base["km1_extension_rows"] != len(covers) * NC: fail("extension row total mismatch")
log(f"BASE FAMILY {json.dumps(base)}")
log(f"TIMING {json.dumps({k: round(v, 1) for k, v in timing.items()})}")

# ---------- cascade tally + persistent rows ----------
tally = {}; persistent = []; l2_total = 0
for f in glob.glob(f"{OUT}/filt/*.out"):
    for line in open(f):
        if not line.startswith("SURVIVOR"): continue
        tally["survivor_rows"] = tally.get("survivor_rows", 0) + 1
        l2_total += int(re.search(r" l2=(\d+)", line)[1])
        if "PERSISTENT" in line:
            tally["persistent"] = tally.get("persistent", 0) + 1; persistent.append(line.strip())
        else:
            for lv, cnt in re.findall(r"l(\d+)=(\d+)", line)[1:]:
                if cnt == "0":
                    tally[f"died_l{lv}"] = tally.get(f"died_l{lv}", 0) + 1; break
tally["level2_improper_lifts_total"] = l2_total
log(f"CASCADE {json.dumps(tally)}")

FOLD = [min(x, P - x) for x in range(P)]
def canon(cls):
    best = None
    for g in range(1, P):
        img = tuple(sorted(FOLD[g * c % P] for c in cls))
        if best is None or img < best: best = img
    return best

ROW_RE = re.compile(r"SURVIVOR base:((?: \d+){%d})" % K)
def canon_row(line):
    return canon([int(x) for x in ROW_RE.match(line).group(1).split()])

# Small gates (general variant) leave millions of persistent rows; the pure-Python
# canonicalisation is the only single-threaded stage, so it is spread over a fork
# pool (Linux). Same function, same result, deterministic order.
t0 = time.time()
if len(persistent) > 20000 and hasattr(os, "fork") and WORKERS > 1:
    import multiprocessing
    with multiprocessing.get_context("fork").Pool(WORKERS) as pool:
        keys = pool.map(canon_row, persistent, chunksize=4096)
else:
    keys = [canon_row(line) for line in persistent]
orbits = {}                      # canonical representative -> number of persistent rows in the orbit
for key in keys:
    orbits[key] = orbits.get(key, 0) + 1
n_persistent = len(persistent)
del keys, persistent             # millions of rows: free them before the kills (memory)
timing["orbit_canon_wall_secs"] = round(time.time() - t0, 1)
reps_sorted = sorted(orbits)   # the complete representative list is written to SUMMARY.json; the log shows the first 20
log(f"PERSISTENT rows={n_persistent} orbits={len(orbits)} (canon {timing['orbit_canon_wall_secs']}s): "
    + "; ".join(" ".join(map(str, k)) for k in reps_sorted[:20]) + (" ; ..." if len(reps_sorted) > 20 else ""))

if shard_n > 1:
    json.dump(dict(p=P, shard=[shard_i, shard_n], precondition=pre, base_family=base, cascade=tally,
                   persistent_orbits=[" ".join(map(str, k)) for k in sorted(orbits)], timing=timing,
                   ir_jobs=ir_res, km1_jobs=km1_res, ext_jobs=ext_res), open(f"{OUT}/SHARD_SUMMARY.json", "w"), indent=1)
    log("SHARD_DONE"); sys.exit(0)

# ---------- 6. kills ----------
kills = []
if not NO_KILL:
    reps = sorted(orbits)
    KILL_RE = r"witness_free_level_(\d+)=(-?\d+) improper_after_gcd=(-?\d+) improper_after_neargcd=(-?\d+)"
    def kill_job(idx):
        # full mode first (exact counts of all witness-free lifts, gcd and near-gcd
        # verdicts); if it times out or overflows the solution cap (small p: the
        # settled CRT family is huge) rerun in prune mode, which drops settled
        # lifts during the search and reports the unsettled ones exactly.
        # Resumable: a kill file with a complete verdict line is reused (a restarted
        # gate with millions of orbits must not redo finished kills). Only the parsed
        # verdict is kept in memory; the full tool output stays in the kill file.
        rep = reps[idx]; kf = f"{OUT}/kills/kill_orbit{idx}.txt"
        if os.path.exists(kf):
            txt = open(kf).read()
            if re.search(KILL_RE, txt):
                return parse_kill(idx, rep, txt, "", 0, 0.0, "prune" if txt.startswith("[full mode failed") else "full")
        cmd = [KILL, str(P)] + [str(c) for c in rep]
        out, err, rc, secs = run(cmd, timeout=KILL_TIMEOUT)
        mode = "full"
        if not re.search(KILL_RE, out):
            out2, err2, rc2, secs2 = run([KILL, "prune", str(P)] + [str(c) for c in rep], timeout=KILL_TIMEOUT)
            out = f"[full mode failed: rc={rc} {(err or '')[-120:].strip()}]\n" + out2; err = err2; rc = rc2; secs += secs2; mode = "prune"
        open(kf, "w").write(out + err)
        return parse_kill(idx, rep, out, err, rc, secs, mode)
    def parse_kill(idx, rep, out, err, rc, secs, mode):
        m = re.search(KILL_RE, out)
        if not m:   # timeout, solution cap or crash in both modes: counted as not killed
            k = dict(orbit=" ".join(map(str, rep)), level=K + 1, witness_free=-1, improper_after_gcd=-1, improper_after_neargcd=-1,
                     mode=mode, secs=round(secs, 2), rows=orbits[rep], error=(err or "")[-200:].strip(), rc=rc)
            return idx, rep, k, rc, False
        k = dict(orbit=" ".join(map(str, rep)), level=int(m[1]), witness_free=int(m[2]), improper_after_gcd=int(m[3]),
                 improper_after_neargcd=int(m[4]), mode=mode, secs=round(secs, 2), rows=orbits[rep])
        ms = re.search(r"witness_free_level_(\d+)=(\d+) witness_free_level_(\d+)=(\d+) crt_candidates=(\d+)", out)
        if ms: k.update({f"witness_free_level_{ms[1]}": int(ms[2]), f"witness_free_level_{ms[3]}": int(ms[4]), "crt_candidates": int(ms[5])})
        return idx, rep, k, rc, True
    t_kill = time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:      # kills are independent single-thread processes
        results = list(ex.map(kill_job, range(len(reps))))
    bad = []
    for idx, rep, k, rc, ok in results:
        kills.append(k)
        if not ok:
            bad.append(rep); log(f"KILL orbit {idx} FAILED rc={rc} {k['error']}"); continue
        if len(reps) <= 20 or k["improper_after_neargcd"] != 0 or rc != 0: log(f"KILL orbit {idx} {k}")
        if k["improper_after_neargcd"] != 0 or rc != 0: bad.append(rep)
    timing["kill_wall_secs"] = round(time.time() - t_kill, 1)
    n_gcd_only = sum(1 for k in kills if k["improper_after_gcd"] == 0)
    n_neargcd = sum(1 for k in kills if k["improper_after_gcd"] != 0 and k["improper_after_neargcd"] == 0)
    log(f"KILLS {len(reps)} orbits: {n_gcd_only} closed by grid witnesses + gcd clause, {n_neargcd} need the near-gcd shift lemma, "
        f"{len(bad)} not killed; cpu={sum(k['secs'] for k in kills):.1f}s wall={timing['kill_wall_secs']}s")
    if bad:
        fail(f"IMPROPER LIFT SURVIVES at level {K+1} after gcd and near-gcd for {len(bad)} orbit(s), first {bad[0]} - GATE NOT CLOSED")
# K=15: K+1=16 is a power of two, so the binary cascade (levels 2..LEVEL) already
# contains the final level 16; a gate closes iff no row persists to level LEVEL.
closed = ("covers_upto_smax" in pre) and (pre["covers_upto_smax"] == 0 or GENERAL) and len(orbits) == 0
NEARGCD_USED = any(k["improper_after_gcd"] != 0 for k in kills)   # some orbit is settled only analytically
SHIFT_USED = base["lifts_settled_by_shift_lemma"] > 0

# ---------- 7. summary + manifest ----------
manifest = {}
for pat in [f"{OUT}/rows_ir/*.txt", f"{OUT}/rows_ir/*.log", f"{OUT}/rows_km1/*.txt", f"{OUT}/rows_km1/*.log",
            f"{OUT}/filt/*.out", f"{OUT}/filt/*.stats", f"{OUT}/kills/*.txt", f"{OUT}/km1low_{P}.*", ENGINE, CASCADE, KILL]:
    for f in sorted(glob.glob(pat)):
        manifest[os.path.relpath(f, OUT).replace("\\", "/")] = sha256(f)
HEAD = (f"Assuming LRC(m) for all m < {K} (arXiv:2609.02604v2), every primitive counterexample v to LRC({K}) "
        f"({K} speeds, {K+1} runners) with pairwise distinct coordinates has {P} | v_1...v_{K}. ")
if SHIFT_USED:
    claim = (HEAD +
             f"Every tuple in I({K},{P},1) has, at level {LEVEL}, only lifts that are ({K},{P},{LEVEL})-proper (witness on the "
             f"1/({LEVEL}*{P}) grid or the gcd clause) or satisfy the power-of-two shift lemma (some d = 2^j, 2 <= j <= {JMAX}, "
             f"with sum over coordinates not divisible by d of g*ceil(2(d/g)/{K+1}) < d, g = gcd(d, coordinate); uses LRC "
             f"for at most {K-2} speeds). {base['lifts_settled_by_shift_lemma']} lifts use the lemma, so this is the "
             f"divisibility of Lemma 2.2(ii) of arXiv:2609.02604v2, not a J({K},{P}) = {{}} certificate.")
    claim_type = "divisibility (framework + power-of-two shift lemma)"
elif NEARGCD_USED:
    claim = (HEAD + f"Some orbits use the near-gcd shift lemma at level {K+1} (K=14 pipeline only).")
    claim_type = "divisibility (framework + near-gcd shift lemma)"
else:
    claim = (HEAD + f"Every tuple in I({K},{P},1) is eventually ({K},{P})-proper; hence J({K},{P}) is empty "
             f"(Sungkawichai-Trakulthongchai, Def 2.1/2.3, Lemma 2.4, Prop 5.1), shown by the binary cascade to level {LEVEL} "
             f"(K+1 = {K+1} is itself a binary level), and Lemma 2.2(i) of arXiv:2609.02604v2 gives the divisibility.")
    claim_type = "J(K,p) empty (framework only)"
summary = dict(
    gate=dict(k=K, p=P, lift_level=LEVEL),
    claim=claim, claim_type=claim_type,
    status="GATE_CLOSED" if closed else "NOT_CLOSED",
    precondition=pre, base_family=base, cascade=tally,
    persistent_orbits=dict(count=len(orbits), representatives=[" ".join(map(str, k)) for k in sorted(orbits)],
                           rows_per_orbit={" ".join(map(str, k)): v for k, v in orbits.items()}),
    level_kills=kills, timing={k: round(v, 1) for k, v in timing.items()},
    tools=dict(engine=os.path.basename(ENGINE), cascade=os.path.basename(CASCADE), kill=os.path.basename(KILL)),
)
json.dump(summary, open(f"{OUT}/SUMMARY.json", "w"), indent=2)
json.dump(manifest, open(f"{OUT}/MANIFEST_SHA256.json", "w"), indent=2)
log(f"{summary['status']} p={P}: manifest entries={len(manifest)}")
if closed: open(f"{OUT}/GATE_CLOSED", "w").write(time.strftime("%Y-%m-%dT%H:%M:%S\n"))
sys.exit(0 if closed else 2)
