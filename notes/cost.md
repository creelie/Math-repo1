# Measured cost

All times are from the 4-core container (one CPU per job unless stated).
"ir gen" is the generation of irredundant 15-class covers, "km1" the
reducible branch (14-class covers extended by every class), and the cascade
times are those of the binary lift cascade.

| p | variant | covers on ≤ 13 classes | ir gen (CPU) | km1 gen (CPU) | cascade | total |
| --- | --- | --- | --- | --- | --- | --- |
| 31 | decomposition | 0 | under 1 s | under 1 s | 35 s | 35 s |
| 239 | decomposition | 0 | 6,921 s | 337 s | 22 s | about 2.0 CPU-hours |
| 179 | general | 83,239 | 4,448 s | 765 s (2,676,385 covers) | 5,441 s irredundant, 24,952 s reducible (summed job time on a shared machine) | about 10 hours of summed job time |
| 269 | decomposition | 0 | 17,222 s | 877 s | 166 s | about 5.1 hours of summed job time (two workers on a shared machine) |
| 271 | decomposition | 0 | 14,093 s | 565 s | 97 s | about 4.1 hours of summed job time (two workers on a shared machine) |
| 281 | decomposition | 0 | 25,828 s | 1,112 s | 90 s | about 7.5 hours of summed job time (two workers on a shared machine) |
| 283 | decomposition | 0 | 23,545 s | 940 s | 68 s | about 6.9 hours of summed job time (two workers on a shared machine) |

The classification run (`km1low 13`) counts the covers on at most 13
classes. It found none for p = 269, 271, 281, 283, 293, 307, 311, 313, 317
and 331, so these use the cheap decomposition variant (163 to 526 s per
prime above 271). There are 26 such covers at p = 251, 39 at 263, 65 at 277,
715 at 257 and 936 at 241.

The generation time grows with p roughly like the number of irredundant
covers; at K = 15 the decomposition primes cost about five times as much as
at K = 14.

Small primes vary. The whole tree of a p = 23 row, to level 128, takes
about ten seconds in the cascade and about a hundred seconds in the
independent brute-force check `check/lift_check_sym.cpp`, and the 91
orbits of p = 23 took 56 minutes on three cores. At p = 19 a single row
can take hours (see `small_primes.md`), so 19 and 17 cost far more per unit
of log-sum than the large primes of the decomposition variant.
