# Measured cost

All times are from the 4-core container (one CPU per job unless stated).
"ir gen" is the generation of irredundant 15-class covers, "km1" the
reducible branch (14-class covers extended by every class), and the cascade
times are those of the binary lift cascade.

| p | variant | covers on ≤ 13 classes | ir gen (CPU) | km1 gen (CPU) | cascade | total |
| --- | --- | --- | --- | --- | --- | --- |
| 31 | decomposition | 0 | under 1 s | under 1 s | 35 s | 35 s |
| 239 | decomposition | 0 | 6,921 s | 337 s | 22 s | about 2.0 CPU-hours |
| 179 | general | 83,239 | 4,464 s | 765 s (2,676,385 covers) | about 20 CPU-minutes per extension shard, 8 shards | about 4.5 CPU-hours |

The classification run (`km1low 13`, at most 150 s per prime) found no
cover on at most 13 classes for p = 269, 271 and 283, so these use the
cheap decomposition variant; p = 251 has 26 such covers and p = 263 has 39.
The primes 257, 277, 281, 293, 307 and 311 did not finish the
classification in 150 s.

The generation time grows with p roughly like the number of irredundant
covers; at K = 15 the decomposition primes cost about five times as much as
at K = 14.

Small primes are cheap: the whole tree of a p = 23 row, to level 128, takes
about ten seconds in the cascade and about a hundred seconds in the
independent brute-force check `check/lift_check_sym.cpp`.
