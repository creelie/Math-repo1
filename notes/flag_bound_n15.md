# The flag bound for sixteen runners

The 16-runner conjecture LRC(15) concerns 15 speeds. Allikvere's flag bound
(arXiv:2609.02604v2, Theorem 3.8) is stated for every n, with two
inputs that depend on n: the ratio hypothesis c_{i+1}/c_i > 4/3, and a shape
bound R_n < r_n on the positive simplex (Lemma 3.7 there, proved for n = 13
and n = 14). This note supplies both for n = 15.

## Statement

**Theorem.** Assume LRC(m) for all m < 15 (Allikvere). If v ∈ Z^15_{>0} is a
primitive counterexample to LRC(15), then

    log(v_1 ⋯ v_15) < 15 log(A_15 / 30) = 497.0298125704…,

with A_15 = K_15^{-1/2} and

    K_15 = 233570089 / 12695029890685690413859209216000000000000.

Every such counterexample has distinct coordinates (LRC(14) handles repeated
speeds), and 720720 = lcm(2, …, 16) divides v_1 ⋯ v_15 (Lemma 3.10 there). So
if P is a set of primes p ≥ 17 with p | v_1 ⋯ v_15 for every primitive
counterexample with distinct coordinates, and the sum of log p over P exceeds
483.5418065791…, then LRC(15) holds.

## Ratio hypothesis

With B_r = (2r/(16(16 − r)))^2 and c_r = B_r − B_{r−1} for 1 ≤ r ≤ 14, the
least ratio c_{r+1}/c_r is 3751/2349 > 4/3. `check/bound_n15.py` computes
K_15 = ∏ c_r and the constants above in exact arithmetic.

## Shape bound: R_15 < 1/2

**Lemma.** On the positive simplex, R_15 < 1/2.

*Proof.* We follow the proof of Lemma 3.7 of Allikvere at n = 15. Put

    q(u) = 1 + 2u − u²,   a(u) = q(u)/u^{2/15},   w(u) = u²/q(u),
    H(u) = (Σ_i w(u_i)) ∏_i a(u_i),

where u_i = x_i / max_j x_j. Then R_15(x)² = 225/H(u), so it suffices to show
H > 900 on (0, 1]^15 with some u_i = 1. Since w(1) = 1/2, w ≥ 0 and
a(u) → ∞ as u → 0, the function H tends to infinity when a coordinate tends
to 0, so it attains its minimum.

*Interior coordinates.* Fix all coordinates but one, u ∈ (0, 1), and put
C = Σ_{j≠i} w(u_j) ∈ (0, 7]. Up to the positive factor ∏_{j≠i} a(u_j), the
objective is φ(u) = a(u)(C + w(u)), and

    φ'(u) = −(2a(u) / (15 u q(u))) β(u),   β(u) = C g(u) − 14u²,
    g(u) = 14u² − 13u + 1.

Let α = (13 − √113)/28 ∈ (0.0846, 0.0847) be the smaller root of g. Then
β(0) = C > 0 and β(α) = −14α² < 0. If C ≤ 1, β is concave with β(0) > 0 and
has one positive root. If C > 1, β is convex with β(1) = 2C − 14 ≤ 0, so its
second root is at least 1. In both cases β has exactly one root u_0 in (0, 1),
it lies in (0, α), and φ decreases on (0, u_0) and increases on (u_0, 1).
Hence every interior coordinate of a minimiser is a root of β in (0, α). With
A = C + w(u) = Σ_i w(u_i), the equation β(u) = 0 reads

    L(u) := q(u) g(u) / (15 u² (1 + u)) = 1/A,

and

    L'(u) = −2P(u) / (15 u³ (1 + u)²),   P(u) = 7u⁵ + 14u⁴ − 27u³ − 11u² − 4u + 1.

On (0, α), P(u) > 1 − 4u − 11u² − 27u³ > 1 − 4(0.0847) − 11(0.0847)² −
27(0.0847)³ > 0, so L is strictly decreasing there. Since A is the same for
all coordinates, all interior coordinates of a minimiser share one value t.

*Two or more coordinates equal to one.* The logarithmic derivative of a is
−2g(u)/(15u q(u)), so a decreases on (0, α), increases up to the larger root
of g and then decreases to a(1) = 2. Hence a ≥ min(a(α), 2) on (0, 1], and

    a(α) > q(0.0846) / 0.0847^{2/15} > 8/5,

an exact comparison after raising to the fifteenth power. If k ≥ 2
coordinates equal one, then Σ w(u_i) ≥ k/2 and
H ≥ (k/2) 2^k (8/5)^{15−k}. This increases with k, and at k = 2 it is
4 (8/5)^{13} > 1801 > 900.

*Exactly one coordinate equal to one.* The other fourteen equal t ∈ (0, α), and

    H(t) = (1/2 + 14t²/q(t)) · 2 · q(t)^{14} / t^{28/15}
         = (1 + 2t + 27t²) q(t)^{13} / t^{28/15}.

Its logarithmic derivative is

    28 (t − 1)(378t³ + 25t² + 10t − 1) / (15 t (t² − 2t − 1)(27t² + 2t + 1)),

which on (0, 1) has the sign of the increasing cubic 378t³ + 25t² + 10t − 1.
So H decreases up to the root t_0 ∈ (0.0724, 0.0725) of the cubic and
increases afterwards, and t_0 < α. Since q and 1 + 2t + 27t² increase and
t^{−28/15} decreases,

    H(t_0) ≥ (1 + 2(0.0724) + 27(0.0724)²) q(0.0724)^{13} / 0.0725^{28/15} > 942.5 > 900,

an exact comparison after raising to the fifteenth power. So H > 900 in every
case, and R_15 < 15/30 = 1/2. ∎

`check/shape_n15.py` checks every identity above symbolically and every
numerical comparison in exact rational arithmetic. The true maximum of R_15 is
about 0.48796, so r_15 = 1/2 loses about 0.37 in the logarithmic bound; we use
1/2 for a clean statement.
