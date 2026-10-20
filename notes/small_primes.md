# Small primes

The gates use primes p ≥ 17 (the primes 2, …, 13 enter through the forced
divisor 720720). This note records how the small ones close.

## The depth of the lift cascade

Allikvere's cascade lifts a level-one cover to the binary levels 2, 4, …, 32
and stops there. At K = 14 the next step is the level-15 fibre split, so
level 32 is enough. At K = 15 the final level 16 is itself binary, and the
port first kept the depth 32. At the small primes this cap was the only
obstacle: rows in which five to seven speeds share a class keep improper
lifts at level 32, but every such row tested dies at level 64 or 128.

| p | row (speed classes) | improper lifts at 16, 32, 64, 128 |
| --- | --- | --- |
| 19 | 1×7, 2, …, 9 | 3,984,537,312; 11,836,500,480; 21,372,422,400; 0 |
| 23 | 1×5, 2, …, 11 | 21,976,320; 593,280; 0 |
| 23 | 1×3, 2×2, 3×2, 4, …, 11 | 25,459,728; 3,335,616; 50,688; 0 |
| 37 | 1×7, 3, 4, 7, 9, 10, 11, 12, 16 | 29,756,160; 82,615,680; 0 |
| 41 | 1×6, 3, 4, 7, 10, 11, 12, 13, 16, 18 | 7,687,200; 11,865,600; 0 |

The counts are labelled lifts (every coordinate distinguishable). Going
deeper is sound for the same reason as the cascade itself: a lift proper at
level L stays proper at level 2L (a witness T/(Lp) is 2T/(2Lp)), the gcd
clause depends only on parities, and the shift lemma at level L only uses
residues modulo d | L. The cascade now runs to level 128 by default
(`MAX_LIFT_EXP_DEF`, printed as `maxlevel`) and the gate's claim names that
level. The independent checker `check/lift_check_sym.cpp` recomputes the
rows above from the definitions, with every time 1, …, Lp − 1 tested and its
own symmetry reduction, and gives the same counts at every level and the
same number of lifts settled by the shift lemma.

## Primes below 32: every class is needed

For p < 32 a speed class c covers only the times T ≡ ±c⁻¹ (mod p), since
16 · dist(Tc, p) < p forces dist = 1. A level-one cover therefore contains
every one of the (p − 1)/2 classes, and the remaining 15 − (p − 1)/2 speeds
are free. `src/small_gate.py` checks this on the covering relation, lists
all such multisets, keeps one representative of each orbit under
multiplication by units modulo p (multiplication by a unit U ≡ 1 modulo the
binary level maps the lift tree of a row onto that of its image and
preserves witnesses, parities and 2-adic valuations), and runs the cascade.
It does not use the cover generator. The representatives are run in
chunks of at most four rows; a finished chunk is kept, so an interrupted
gate resumes where it stopped. The gates at 23, 29 and 31 were run before
the chunking was added (commit 2bdce84, one file per worker); their
manifests record that version of the script.

| p | classes | free speeds | multisets | orbits | result |
| --- | --- | --- | --- | --- | --- |
| 31 | 15 | 0 | 1 | 1 | closed, every lift dies by level 16 |
| 29 | 14 | 1 | 14 | 1 | closed, every lift dies by level 32 |
| 23 | 11 | 4 | 1,001 | 91 | closed; 57 orbits die at level 64 and 34 at level 128, none earlier; 8,636,800 lifts by the shift lemma |
| 19 | 9 | 6 | 3,003 | 335 | not run to the end: too costly here (see below) |
| 17 | 8 | 7 | 3,432 | 429 | not run: more free speeds than 19 |

At p = 19 single rows are far heavier than at 23. Measured on one core with
the cascade capped at lower levels: the row 1⁷,2,…,9 reaches level 16 in
2.8 s, level 32 in 16 s and level 128 in 284 s; the row 1⁴,2,3³,4,5,6²,7,8,9
needs 89 s for level 16 alone and more than 300 s for level 32; the row
1³,2³,3,4³,5,6,7,8,9 needs 34 s for level 16. Two chunks of four such rows
ran 70 CPU-minutes each without finishing. With the cost growing about
fourfold per doubling of the level, one row can take hours, and the 335
orbits of p = 19 would cost of the order of a hundred CPU-hours for
log 19 = 2.94 — far more per unit of log-sum than a large prime of the
decomposition variant (about two to four CPU-hours for log p ≈ 5.6). The
gates at 19 and 17 are therefore left for later; a run there needs the
lift tree of a row split into subtrees, since one row can outlast a
two-hour job.

## Primes from 37 on

From p = 37 on a cover may omit classes, so these primes go through the full
gate (`src/gate_k15.py`): the reducible branch is the general `km1all`
family of 14-class covers, each extended by every class.
