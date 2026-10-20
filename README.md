# Sixteen lonely runners: prime gates for fifteen speeds

This repository carries the modular method that proved the Lonely Runner
Conjecture for fifteen runners (J. Allikvere, *Fourteen and fifteen lonely
runners*, arXiv:2609.02604v2) one step further, to sixteen runners, and
records which primes it settles there.

The conjecture LRC(k) says that for any k distinct positive integer speeds
v_1, ..., v_k there is a time t with ||t v_i|| ≥ 1/(k+1) for every i, where
||x|| is the distance to the nearest integer. Sixteen runners is LRC(15). It
is open.

## What is proved here

Assume LRC(m) for all m ≤ 14 (Allikvere's theorem). For every prime p in the
table below, every primitive counterexample v to LRC(15) with pairwise
distinct coordinates has

    p | v_1 v_2 ... v_15.

Together with Rosenfeld's forced divisor lcm(2, ..., 16) = 720720, the product
of the speeds of a counterexample is divisible by 720720 times the product of
these primes.

| p | gate | lift level | rows through the cascade | lifts settled by the shift lemma | summary |
| --- | --- | --- | --- | --- | --- |
| 23 | small-prime | 128 | 91 orbits (1,001 covers) | 8,636,800 | `results/summary_23.json` |
| 29 | small-prime | 128 | 1 orbit (14 covers) | 0 | `results/summary_29.json` |
| 31 | small-prime | 128 | 1 orbit (1 cover) | 0 | `results/summary_31.json` |
| 179 | general | 32 | 438,266,496 rows | 403,094,340 | `results/summary_179.json` |
| 239 | decomposition | 32 | 10,536,820 rows | not used | `results/summary_239.json` |
| 269 | decomposition | 128 | 6,596,721 rows | not used | `results/summary_269.json` |
| 271 | decomposition | 128 | 4,179,959 rows | not used | `results/summary_271.json` |

This does not prove LRC(15). Allikvere's flag bound, at n = 15, gives
log(v_1 ... v_15) < 497.03 for a primitive counterexample, so the method
closes the conjecture once the verified primes p ≥ 17 have
Σ log p > 483.54, for instance every prime from 17 to 523 (93 primes). The
bound needs a shape estimate R_15 < 1/2 that the paper proves only for
n = 13 and 14; it is proved here for n = 15 (`notes/flag_bound_n15.md`,
checked exactly by `check/shape_n15.py`), and `check/bound_n15.py`
recomputes both constants. The primes settled so far are far from that sum.

## What is new

1. **The port to fifteen speeds.** Allikvere's last step works at level
   15 = 3 · 5. For fifteen speeds that level is 16 = 2^4, a binary level, so
   the factorisation is unavailable. Allikvere lists this as future work.
2. **A power-of-two shift lemma** (`notes/shift_lemma.md`). For d = 2^j, if
   the speeds not divisible by d satisfy
   Σ gcd(d, v_i) · ⌈2 (d / gcd(d, v_i)) / 16⌉ < d, the tuple has a lonely
   time. It settles lifts that no witness time settles.
3. **Deeper lifts** (`notes/small_primes.md`). The cascade now lifts to
   level 128 instead of 32. At the small primes, rows with five to seven
   speeds in one class keep improper lifts at level 32 and all die at level
   64 or 128.
4. **A gate for the primes below 32** (`src/small_gate.py`). There every
   class is needed in a cover, so the covers are listed directly, without
   the cover generator.
5. **An exact symmetry reduction.** Equal coordinates of a row are
   interchangeable, so only canonical lift patterns are expanded, each
   weighted by a product of binomial coefficients.
6. **The shape bound at n = 15** (`notes/flag_bound_n15.md`): R_15 < 1/2
   on the positive simplex, which completes the flag bound for sixteen
   runners.
7. **Independent checks.** `check/lift_check_sym.cpp` recomputes the lift
   tree of a row to any binary level from the definitions (every time
   tested, its own symmetry reduction, the shift lemma by brute force) and
   agrees exactly with the cascade; `check/compare_cascade.py` compares the
   two row by row. `check/lift_check.cpp` and `check/lift_check_shift.cpp`
   recount small cases without symmetry, and `check/shift_lemma_check.py`
   verifies the lemma's mechanism with exact rationals.

## What limits the method

Nothing found so far blocks a prime outright; the limit is cost. From
p = 37 on the reducible covers form the general `km1all` family (19,324
fourteen-class covers at p = 37, each extended by all 18 classes), and rows
with many repeated classes need the deep lifts, so these primes cost more
than a container can give. Above them every prime needs its irredundant
covers generated, which takes hours per prime (`notes/cost.md`).

## Build and run

    make             # builds bin/ (g++ with C++17 and C++20)
    make selftest    # cascade self-test, shift-lemma check, bound checks
    scripts/small15.sh 2 29 31       # gates for the primes below 32
    scripts/sweep15.sh 4 239 251     # gates for p = 239 and 251 on 4 workers
    python3 check/compare_cascade.py 23 rows.txt   # cascade against the checker

Each gate writes `SUMMARY.json` and `MANIFEST_SHA256.json` (SHA-256 of every
row file, cascade output and binary) to `results/`. The bulky row files are
deleted after each gate unless `KEEP=1`; rerunning a gate regenerates them
exactly.

The gate for p = 239 ran with an earlier build of the cascade, before the
shift lemma, the symmetry reduction and the deeper lifts were added. Every
row died before level 32 and no lift needed the shift lemma, so its claim
(J(15, 239) is empty) does not depend on the later changes. Its claim text
also names a level-16 kill step; no orbit persisted, so that step never ran.

The gate for p = 179 ran with the level-32 build of the cascade with the
shift lemma, before the deeper lifts were added. Every row died by level 16
(101,275,180 at level 4, 2,679,532 at level 8, 2,045 at level 16), so its
claim does not depend on them. Its claim text writes LRC(16) for sixteen
runners; `claim_note` in the summary says so. At p = 179 there are covers
on 13 classes, so the reducible branch is the general one: 2,676,385
fourteen-class covers, each extended by all 89 classes.

## Computation

All gates ran on a 4-core container. GitHub Actions is not used for the
campaign. Measured costs are in `notes/cost.md`.

## Credits and licence

The cover generator and the cascade are Allikvere's code (CC BY 4.0,
doi:10.5281/zenodo.22667683); `NOTICE.md` lists exactly what was changed.
The framework is due to Sungkawichai and Trakulthongchai (arXiv:2604.23906)
and the forced divisor to Rosenfeld. Everything else is under the MIT
licence.
