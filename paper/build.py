#!/usr/bin/env python3
"""Assemble the two versions of the paper from paper/src and results/.

usage: python3 paper/build.py [DATE]

Writes paper/hexagon/main.tex (REVTeX, two columns) and paper/journal/main.tex
(amsart), each with a copy of paper/figures/*.png. Every number that depends
on the gates (the primes of Theorem 2, the tables, the log-sum) is read from
the summaries results/summary_<p>.json of the closed gates.
"""
import glob, json, math, os, shutil, sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATE = sys.argv[1] if len(sys.argv) > 1 else "19 October 2026"
SRC = os.path.join(HERE, "src")


def thousands(n):
    return "{,}".join(f"{n:,}".split(","))


def num(n):
    return "$" + thousands(n) + "$"


# ------------------------------------------------------------- gate summaries
gates = {}
for f in glob.glob(os.path.join(ROOT, "results", "summary_*.json")):
    s = json.load(open(f))
    if s.get("status") == "GATE_CLOSED" and int(s["gate"]["p"]) >= 17:
        gates[int(s["gate"]["p"])] = s
primes = sorted(gates)
if not primes:
    sys.exit("no closed gate in results/")


def died_level(s):
    lv = [int(k[6:]) for k, v in s["cascade"].items() if k.startswith("died_l") and v]
    return max(lv) if lv else 1


def shift_count(s):
    if "lifts_settled_by_shift_lemma" in s.get("cascade", {}):
        return s["cascade"]["lifts_settled_by_shift_lemma"]
    return s.get("base_family", {}).get("lifts_settled_by_shift_lemma", 0)


def hours(s):
    if "jobs" in s:
        return sum(j["secs"] for j in s["jobs"]) / 3600
    t = s["timing"]
    return (t["ir_gen_cpu_secs"] + t["km1_gen_cpu_secs"] + t["ir_cascade_wall_secs"]
            + t["km1_cascade_wall_secs"] + t.get("precondition_secs", 0)) / 3600


def fmt_hours(h):
    return "$<0.1$" if h < 0.1 else f"${h:.1f}$"


small = [p for p in primes if "covers" in gates[p]]
large = [p for p in primes if p not in small]

small_rows = []
for p in small:
    s = gates[p]
    c = s["covers"]
    small_rows.append(f"${p}$ & ${c['free_speeds']}$ & {num(c['multisets'])} & "
                      f"{num(c['orbit_representatives'])} & level ${died_level(s)}$ & "
                      f"{num(shift_count(s))}\\\\")

result_rows = []
for p in primes:
    s = gates[p]
    if p in small:
        variant, vectors = "small prime", s["covers"]["orbit_representatives"]
    else:
        bf = s["base_family"]
        variant = "decomposition" if bf["reducible_variant"] == "decomposition" else "general"
        vectors = bf["irredundant_rows"] + bf["km1_extension_rows"]
    result_rows.append(f"${p}$ & {variant} & {num(vectors)} & ${died_level(s)}$ & "
                       f"{num(shift_count(s))} & {fmt_hours(hours(s))}\\\\")

logsum = sum(Decimal(p).ln() for p in primes)
product = math.prod(primes)
plist = ",".join(str(p) for p in primes)
words = ", ".join(str(p) for p in primes[:-1]) + (" and " if len(primes) > 1 else "") + str(primes[-1])

# shift-lemma remark
used = [p for p in primes if shift_count(gates[p])]
unused = [p for p in primes if not shift_count(gates[p])]
parts = [f"at $p={p}$ it settles {num(shift_count(gates[p]))} lifts" for p in used]
if parts:
    shiftremark = ("The shift lemma is used in the proof: " + "; ".join(parts))
    if unused:
        shiftremark += (", while at " + ", ".join(f"${p}$" for p in unused[:-1])
                        + (" and " if len(unused) > 1 else "") + f"${unused[-1]}$"
                        + " no lift needs it")
    shiftremark += ". Without it these gates would leave improper lifts that no witness on the binary grid settles."
else:
    shiftremark = ""

# depth remark from the small-prime gates
dr = []
for p in small:
    cs = gates[p]["cascade"]
    deep = {int(k[6:]): v for k, v in cs.items() if k.startswith("died_l") and v}
    if max(deep) > 32:
        n = gates[p]["covers"]["orbit_representatives"]
        bits = [f"{thousands(v)} at level ${l}$" for l, v in sorted(deep.items())]
        allbeyond = min(deep) > 32
        dr.append(f"At $p={p}$ {'every one' if allbeyond else 'some'} of the {thousands(n)} orbit "
                  f"representatives {'keeps' if allbeyond else 'keep'} improper lifts at level $32$; "
                  f"they lose the last ones " + ", ".join(bits[:-1]) + (" and " if len(bits) > 1 else "") + bits[-1] + ".")
depthremark = " ".join(dr)

# constants of Table 5
n = 15
B = [F(0)] + [F(2 * r, (n + 1) * (n + 1 - r)) ** 2 for r in range(1, n)]
c = [B[r] - B[r - 1] for r in range(1, n)]
crrows = []
for r in range(1, n):
    cr = c[r - 1]
    if r < n - 1:
        q = c[r] / c[r - 1]
        ratio = f"${q.numerator}/{q.denominator}$ & ${float(q):.4f}$"
    else:
        ratio = " & "
    crrows.append(f"${r}$ & ${cr.numerator}/{cr.denominator}$ & {ratio}\\\\")

abstract = (
    "We prove that the speeds of any counterexample to the lonely runner "
    "conjecture for sixteen runners with no common factor have a product below "
    "an explicit bound and "
    f"divisible by 720720 and by each of the primes {words}. The conjecture for "
    "sixteen runners therefore follows from a finite list of prime checks of the "
    "same kind. Two tools carry the method of the fifteen-runner case one step "
    "further: a shift lemma for powers of two, and binary lifting to greater depth."
)

subst = {
    "\\PAPERDATE": DATE,
    "\\ABSTRACTTEXT": abstract,
    "\\PZEROLIST": plist,
    "\\PZEROPRODUCT": "\\cdot".join(str(p) for p in primes),
    "\\PZEROLOGSUM": f"${logsum:.4f}$".replace("$", "") if False else f"{logsum:.4f}",
    "\\SMALLROWS": "\n".join(small_rows),
    "\\RESULTROWS": "\n".join(result_rows),
    "\\SHIFTREMARK": shiftremark,
    "\\DEPTHREMARK": depthremark,
    "\\COSTONESEVENNINE": f"{hours(gates[179]):.0f}" if 179 in gates else "ten",
    "\\CRROWS": "\n".join(crrows),
}


def fill(text):
    for k in sorted(subst, key=len, reverse=True):
        text = text.replace(k + "{}", subst[k]).replace(k, subst[k])
    return text


body = open(os.path.join(SRC, "body.tex")).read()
tail = [open(os.path.join(SRC, f)).read() for f in ("declarations.tex", "appendix.tex", "refs.tex")]
for kind in ("hexagon", "journal"):
    head = open(os.path.join(SRC, f"head_{kind}.tex")).read()
    out = os.path.join(HERE, kind)
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    for png in glob.glob(os.path.join(HERE, "figures", "*.png")):
        shutil.copy(png, os.path.join(out, "figures"))
    text = fill(head + body + "\n" + "\n".join(tail))
    left = [k for k in subst if k in text]
    if left:
        sys.exit(f"unfilled placeholders: {left}")
    open(os.path.join(out, "main.tex"), "w").write(text)
    print(f"{kind}/main.tex: primes {plist}, log-sum {logsum:.4f}")
