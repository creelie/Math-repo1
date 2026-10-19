# The power-of-two shift lemma

This note states and proves the lemma that replaces the level-15 step of
Allikvere's fifteen-runner proof (arXiv:2609.02604v2, Lemma 4.6,
"Few-exception shift") when the number of speeds is K = 15. Here K + 1 = 16 = 2^4 has no odd
factor, so the factorisation 15 = 3 · 5 used for fourteen speeds is not
available, and the last binary level 32 must settle everything.

Notation: ||x|| is the distance from x to the nearest integer, and LRC(m)
says that for any m distinct positive integers u_1, ..., u_m there is a real
t with ||t u_i|| ≥ 1/(m+1) for all i.

## Statement

**Lemma.** Let K ≥ 3, let u_1, ..., u_K be positive integers, let j ≥ 1 and
d = 2^j. Put E = {i : d ∤ u_i}, and for i in E let g_i = gcd(d, u_i) and
N_i = d / g_i. Assume |E| ≥ 2 and LRC(m) for every m ≤ K − 2. If

    S_d = Σ_{i ∈ E} g_i · ⌈2 N_i / (K + 1)⌉ < d,

then there is a rational t with ||t u_i|| ≥ 1/(K + 1) for every i.

## Proof

Let B be the set of indices with d | u_i. Then |B| = K − |E| ≤ K − 2, and the
numbers u_i / d (i in B) take m ≤ K − 2 distinct values. By LRC(m) there is a
t_1 with ||t_1 u_i / d|| ≥ 1/(m + 1) ≥ 1/(K + 1) for i in B (if B is empty,
take t_1 = 0). Put t_0 = t_1 / d.

Consider the d times t_q = t_0 + q/d, q = 0, ..., d − 1. For i in B,
t_q u_i = t_0 u_i + q (u_i / d) differs from t_0 u_i by an integer, so
||t_q u_i|| = ||t_0 u_i|| ≥ 1/(K + 1) for every q.

Fix i in E. As q runs over 0, ..., d − 1, the fractional part of q u_i / d
runs over the grid (1/N_i) Z / Z, and each grid point is hit exactly g_i
times. The shift q is bad for u_i when t_0 u_i + q u_i / d lies within
distance less than 1/(K + 1) of an integer, that is, when q u_i / d lies in a
fixed open arc of length 2/(K + 1) < 1. An open arc of length 2/(K + 1)
contains at most ⌈2 N_i / (K + 1)⌉ points of a grid of mesh 1/N_i, because an
open interval of length ℓ contains at most ⌈ℓ⌉ integers. So u_i forbids at
most g_i ⌈2 N_i / (K + 1)⌉ of the d shifts.

All exceptional speeds together forbid at most S_d < d shifts, so some t_q
has ||t_q u_i|| ≥ 1/(K + 1) for every i. The set of admissible t_1 is a
finite union of intervals with rational endpoints, so t_1, and hence t_q, can
be taken rational. ∎

## Remarks

1. For d = 2 the sum is S_2 = |E| (each odd speed has g = 1, N = 2 and
   ⌈4/(K+1)⌉ = 1), so the condition S_2 < 2 says that at most one speed is
   odd. This is the gcd clause of the framework with the prime 2, so the
   lemma extends it to higher powers of 2. For an odd prime d it is
   Allikvere's few-exception shift lemma (there g_i = 1 and N_i = d).

2. In the cascade the lemma is applied to a lift w at binary level 2^L, for
   each j = 2, ..., L. The test depends only on the 2-adic valuations of the
   coordinates capped at j, which are constant on the residue class of w
   modulo 2^L p because 2^j divides 2^L p. Hence every integer vector in that
   class with distinct positive coordinates has the lonely-runner property:
   this is the third alternative of Lemma 2.2(ii) (prime-divisibility
   criterion) in Allikvere's paper. A lift settled at level 2^j stays settled
   at every higher level, since its descendants have the same residues
   modulo 2^j p.

3. When the lemma is called, the gcd clause has already failed, so at least
   two coordinates are odd. Hence |E| ≥ 2 automatically and the speeds
   divisible by d number at most K − 2 = 13. Both this input, LRC(m) for
   m ≤ 13, and the hypothesis LRC(m) for all m < 15 of the criterion are
   theorems of Allikvere, *Fourteen and fifteen lonely runners*,
   arXiv:2609.02604v2.

4. In `src/cascade_filter_k.cpp` the function `shift_proper` evaluates S_d
   with g_i = 2^{v_2(u_i)} and ⌈2N/(K+1)⌉ = (2N + K) div (K + 1). The
   self-test compares it with a brute-force count of forbidden shifts on
   20,000 random cases, and `check/lift_check_shift.cpp` recomputes it
   independently on whole rows.

## Example

At p = 179 the binary cascade without the lemma left rows whose surviving
lifts have thirteen coordinates divisible by 32 and two odd ones. For such a
lift and d = 4, the two odd speeds give S_4 = 2 · ⌈8/16⌉ = 2 < 4, so the lift
is settled at level 4.
