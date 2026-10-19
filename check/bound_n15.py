#!/usr/bin/env python3
"""Constants of the flag bound for n = 15 speeds (sixteen runners).

Recomputes, with exact rationals, the constants of Theorem 3.8 (flag bound)
and Lemma 3.10 (forced divisor) of J. Allikvere, "Fourteen and fifteen lonely
runners", arXiv:2609.02604v2, at n = 15:

  B_r = (2r / ((n+1)(n+1-r)))^2,  c_r = B_r - B_{r-1},  K_n = prod c_r,
  A_n = K_n^(-1/2),  log(v_1...v_n) < n log(A_n r_n / n) with r_n = 1/2,

checks the ratio hypothesis c_{r+1}/c_r > 4/3, and subtracts the forced
divisor lcm(2,...,16) = 720720. The shape factor R_15 < 1/2 is proved in
notes/flag_bound_n15.md and checked exactly by check/shape_n15.py.

With an argument DIR it also adds up log p over the primes p whose
SUMMARY file DIR/summary_<p>.json reports a closed gate.
"""
import glob, json, math, os, sys
from fractions import Fraction as F
from decimal import Decimal, getcontext

getcontext().prec = 60
n = 15
d = n - 1
B = [F(0)] + [F(2 * r, (n + 1) * (n + 1 - r)) ** 2 for r in range(1, d + 1)]
c = [B[r] - B[r - 1] for r in range(1, d + 1)]
ratios = [c[i + 1] / c[i] for i in range(d - 1)]
assert min(ratios) == F(3751, 2349) and min(ratios) > F(4, 3)
Kn = F(1)
for ci in c:
    Kn *= ci
# log bound = n log(A_n / (2n)) = n (-(1/2) log K_n - log(2n))
logK = Decimal(Kn.numerator).ln() - Decimal(Kn.denominator).ln()
bound = n * (-logK / 2 - Decimal(2 * n).ln())
lcm = math.lcm(*range(2, n + 2))
assert lcm == 720720 == 2**4 * 3**2 * 5 * 7 * 11 * 13
need = bound - Decimal(lcm).ln()
print(f"n = {n}: min c ratio = {min(ratios)} > 4/3")
print(f"K_15 = {Kn}")
print(f"log(v_1 ... v_15) < 15 log(A_15/30) = {bound:.12f}")
print(f"forced divisor lcm(2..16) = {lcm}, log = {Decimal(lcm).ln():.12f}")
print(f"log-sum needed from primes p >= 17: {need:.12f}")
assert abs(float(bound) - 497.0298125704946) < 1e-9

if len(sys.argv) > 1:
    total, primes = Decimal(0), []
    for f in sorted(glob.glob(os.path.join(sys.argv[1], "summary_*.json"))):
        s = json.load(open(f))
        p = int(s["gate"]["p"])
        if s.get("status") == "GATE_CLOSED" and p >= 17:
            primes.append(p)
            total += Decimal(p).ln()
    primes.sort()
    print(f"closed primes: {primes}")
    print(f"log-sum of closed primes = {total:.6f} of {need:.6f} needed "
          f"({float(total / need) * 100:.1f}%)")
