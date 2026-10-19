#!/usr/bin/env python3
"""Exact checks for the shape bound R_15 < 1/2 (notes/flag_bound_n15.md).

Notation of arXiv:2609.02604v2, Section 3.5, with n = 15:
  q(u) = 1 + 2u - u^2,  a(u) = q(u) / u^(2/n),  w(u) = u^2 / q(u),
  H(u) = (sum_i w(u_i)) * prod_i a(u_i),   R_n^2 = n^2 / H.
R_15 < 1/2 is H > 4 n^2 = 900 on (0,1]^15 with some u_i = 1.
Every comparison below is made in exact rational or integer arithmetic.
"""
from fractions import Fraction as Fr
import sympy as sp

n = 15
u, t = sp.symbols("u t", positive=True)
q = 1 + 2 * u - u**2
g = (n - 1) * u**2 - (n - 2) * u + 1
ok = True

def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok &= bool(cond)

# alpha = smaller root of g = 14u^2 - 13u + 1; it lies in (0.0846, 0.0847).
gs = lambda x: 14 * x * x - 13 * x + 1
check("alpha in (0.0846, 0.0847)", gs(Fr(846, 10000)) > 0 > gs(Fr(847, 10000)))

# phi'(u) has the sign of -beta(u), beta = C g(u) - (n-1) u^2 (derived symbolically).
C = sp.symbols("C", positive=True)
a = q / u**sp.Rational(2, n); w = u**2 / q
phi_prime = sp.diff(a * (C + w), u)
beta = C * g - (n - 1) * u**2
check("phi' = -(2a/(n u q)) beta", sp.simplify(phi_prime + 2 * a / (n * u * q) * beta) == 0)

# Equation (6): beta = 0 with A = C + w(u) is L(u) = 1/A, L = q g / (n u^2 (1+u)).
A = sp.symbols("A", positive=True)
L = q * g / (n * u**2 * (1 + u))
lhs = sp.simplify((beta.subs(C, A - w)) - (g * A - n * u**2 * (1 + u) / q))
check("beta(C=A-w) = g A - n u^2 (1+u)/q, so beta = 0 iff L(u) = 1/A", lhs == 0)
dL = sp.factor(sp.together(sp.diff(L, u)))
P = 7 * u**5 + 14 * u**4 - 27 * u**3 - 11 * u**2 - 4 * u + 1
check("L' = -2 P(u) / (15 u^3 (1+u)^2)", sp.simplify(dL + 2 * P / (15 * u**3 * (1 + u)**2)) == 0)
x = Fr(847, 10000)
check("P > 0 on (0, alpha): 1 - 4x - 11x^2 - 27x^3 > 0 at x = 0.0847", 1 - 4 * x - 11 * x**2 - 27 * x**3 > 0)

# a(u) > 8/5 on (0,1]: log a' has the sign of -(14u^2 - 13u + 1); local minimum at alpha, a(1) = 2.
la = sp.simplify(sp.diff(sp.log(a), u) + 2 * g / (n * u * q))
check("(log a)' = -2 g / (15 u q)", la == 0)
q0 = 1 + 2 * Fr(846, 10000) - Fr(846, 10000)**2
check("a(alpha) > q(0.0846) / 0.0847^(2/15) > 8/5  (15th powers)", q0**15 / Fr(847, 10000)**2 > Fr(8, 5)**15)

# Two or more coordinates equal to one: H >= (k/2) 2^k (8/5)^(n-k), increasing in k.
f = lambda k: Fr(k, 2) * 2**k * Fr(8, 5)**(n - k)
check("k >= 2 ones: (k/2) 2^k (8/5)^(n-k) increasing and > 900 at k = 2",
      all(f(k + 1) > f(k) for k in range(2, n)) and f(2) > 900)

# Exactly one coordinate equal to one; the other fourteen share t in (0, alpha):
#   H(t) = (1 + 2t + 27t^2) q(t)^13 / t^(28/15).
Ht = (q.subs(u, t) + 2 * (n - 1) * t**2) * q.subs(u, t)**(n - 2) * t**sp.Rational(-2 * (n - 1), n)
check("H(t) = (1 + 2t + 27t^2) q^13 t^(-28/15)",
      sp.simplify(Ht - (1 + 2 * t + 27 * t**2) * q.subs(u, t)**13 * t**sp.Rational(-28, 15)) == 0)
cub = 378 * t**3 + 25 * t**2 + 10 * t - 1
dlog = sp.simplify(sp.diff(sp.log(Ht), t) - 28 * (t - 1) * cub / (15 * t * (t**2 - 2 * t - 1) * (27 * t**2 + 2 * t + 1)))
check("(log H)' = 28 (t-1) cubic / (15 t (t^2-2t-1)(27t^2+2t+1))", dlog == 0)
cf = lambda x: 378 * x**3 + 25 * x**2 + 10 * x - 1
check("cubic root t0 in (0.0724, 0.0725), below alpha", cf(Fr(724, 10000)) < 0 < cf(Fr(725, 10000)))
lo, hi = Fr(724, 10000), Fr(725, 10000)
num = (1 + 2 * lo + 27 * lo**2)**15 * (1 + 2 * lo - lo**2)**195
check("H(t0) > 900  (15th powers: increasing factors at 0.0724, t^(-28/15) at 0.0725)",
      num > 900**15 * hi**28)
val = float((num / hi**28)) ** (1 / 15)
print(f"lower bound for H(t0): {val:.3f}  (900 needed)")
print("ALL PASS: R_15 < 1/2" if ok else "SOME CHECK FAILED")
raise SystemExit(0 if ok else 1)
