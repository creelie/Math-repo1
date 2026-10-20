#!/usr/bin/env python3
"""Randomised check of the power-of-two shift lemma for K = 15 speeds.

For random pairwise distinct positive integer tuples u that satisfy the shift
condition for some d = 2^j (j >= 2) -- the speeds not divisible by d forbid, in
total, fewer than d of the shifts q/d -- the script finds a time t0 at which
the speeds divisible by d are at distance >= 1/(m+1) from the integers (m =
their number; the lemma takes t0 from LRC(m), here it is found by search) and
then checks with exact rational arithmetic that some shifted time t0 + q/d is
a lonely-runner witness (every speed at distance >= 1/16).
"""
import random, sys
from fractions import Fraction as F
from math import gcd

K = 15

def dist(x):            # distance of a Fraction to the nearest integer
    x = x - (x.numerator // x.denominator)
    return min(x, 1 - x)

def forbidden(u, d):    # g * ceil(2 (d/g) / (K+1)) for a speed not divisible by d
    g = gcd(u, d); n = d // g
    return g * (-(-2 * n // (K + 1)))

def shift_condition(us, d):
    return sum(forbidden(u, d) for u in us if u % d) < d

def find_t0(D):
    m = len(D)
    if m == 0: return F(0)
    need = F(1, m + 1)
    cands = set()
    for u in D:
        for a in range((m + 1) * u):
            cands.add(F(a, (m + 1) * u))
    for t in sorted(cands):
        if all(dist(t * u) >= need for u in D): return t
    return None

def main(trials, seed):
    rng = random.Random(seed)
    checked = 0
    while checked < trials:
        j = rng.choice([2, 3, 4, 5]); d = 1 << j
        e = rng.randint(1, 7)
        ex = []
        while len(ex) < e:
            u = rng.randint(1, 400)
            if u % d and u not in ex: ex.append(u)
        D = []
        while len(D) < K - e:
            u = d * rng.randint(1, 400 // d + 3)
            if u not in D: D.append(u)
        us = ex + D
        if not shift_condition(us, d): continue
        t0 = find_t0(D)
        assert t0 is not None, ("LRC fails?", D)
        good = [q for q in range(d) if all(dist((t0 + F(q, d)) * u) >= F(1, K + 1) for u in us)]
        bad_shifts = [q for q in range(d) if q not in good]
        assert good, ("no good shift", us, d, t0)
        assert len(bad_shifts) <= sum(forbidden(u, d) for u in ex), ("bound violated", us, d)
        checked += 1
    print(f"shift lemma: {trials} random tuples satisfying the condition, every one has a good shift")

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300, int(sys.argv[2]) if len(sys.argv) > 2 else 1)
