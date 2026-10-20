# Third-party code and credits

## Allikvere's fifteen-runner pipeline (CC BY 4.0)

The level-one cover generator and the binary lift cascade are taken from the
verification record of

> Jaan Allikvere, *Fourteen and fifteen lonely runners*, arXiv:2609.02604v2
> (2026); record *Fifteen lonely runners: manuscript, gate certificates, and
> audit code*, Zenodo, doi:10.5281/zenodo.22667683, licence CC BY 4.0
> (ORCID 0009-0003-5228-7015).

| File here | Upstream file (SHA-256 prefix) | Changes |
| --- | --- | --- |
| `src/basegen_k_campaign.cpp` | `basegen_k_campaign.cpp` (4b48224e9f31550d) | none; built with `-DK=15` |
| `src/bgk_incremental_gain.h` | `bgk_incremental_gain.h` (bd150a98a7fdcf99) | none |
| `src/cascade_filter_k.cpp` | `cascade_filter_k.cpp` (50fa9385abe21075) | `K = 15`; power-of-two shift lemma (`shift_proper`); exact symmetry reduction for equal coordinates with 128-bit weighted counts; lift depth 128 instead of 32 (build option `MAX_LIFT_EXP_DEF`), printed as `maxlevel`; primes 17 to 29 accepted; self-test parts C to E |
| `src/gate_k15.py` | `gate_k14.py` (290f410106541525) | `K = 15`; no level-15 kill (16 is a binary level); shift-lemma accounting; the claim names the lift level reported by the cascade |

These files remain under CC BY 4.0. The changes listed above, and every
other file in this repository, are under the MIT licence in `LICENSE`.

## Framework

The modular framework (proper and improper lifts, eventual properness, the
prime-divisibility criterion) is that of T. Sungkawichai and
T. Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*,
arXiv:2604.23906v2 (2026), in the form of Allikvere's Section 2. The forced
divisor is Rosenfeld's observation (*The lonely runner conjecture holds for
eight runners*, Math. Comp. (2026), doi:10.1090/mcom/4243, Lemma 4).

## Independent checks

`check/lift_check.cpp`, `check/lift_check_shift.cpp`,
`check/lift_check_sym.cpp`, `check/shift_lemma_check.py`,
`check/bound_n15.py` and `check/shape_n15.py` were written for this
repository from the definitions in the papers, without reusing the upstream
sources.
