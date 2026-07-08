# P15-S5 Odd Structure Certificate

Project: Signed Reversal Rank Theorem (ID P15)
Gate: P15-S5 odd-n structure
Date: 2026-07-06
Schema: `p15.s5.odd_structure_certificate.v1`
Status: PASS

## Command

```powershell
python -B scripts\p15_s5_odd_structure_certificate.py --write-json certified\P15_S5_ODD_STRUCTURE_CERTIFICATE.json --limit-m 10 --hm-ground-m 7 --hm-positive-m 16
```

## Checks

| Check | Status | Meaning |
|---|---|---|
| `odd_channel_exact_replay` | PASS | Exact full enumeration for odd `n=3,5,7`; verifies ranks, collapse, `a_c` formula, and trace identity. |
| `odd_center_scalar_sweep` | PASS | Exact even-count sweep for `a_c=d_{k-1}-d_k` through odd `n=21`; empty top rows are fenced. |
| `odd_k1_H_reduction` | PASS | Checks `a(2m+1,1)=a(2m,1)+H_m`; ground residual `m=2..7`; `H_m>0` through `m=16`. |
| `odd_pair_standard_reduction` | PASS | Reduces odd pair-standard channel to exact `lam_prime` and `mu`; verifies eigenvectors for `n=3,5,7`. |

## Key Outputs

Odd exact channel replay includes the base case:

```text
n=3, k=1: rank_M_ref=2, rank_M_pair_std=1, a=-1, a_c=2
```

Thus `n=3,k=1` is rank-correct by nonzero scalar, not by positivity.

For `H_m`, the ground replay confirms:

```text
m=2: a_odd=30, a_even=8, H_m=22
m=3: a_odd=3316, a_even=376, H_m=2940
m=4: a_odd=796744, a_even=62448, H_m=734296
m=5: a_odd=297858864, a_even=17724736, H_m=280134128
m=6: a_odd=161179803488, a_even=7736076800, H_m=153443726688
m=7: a_odd=119371368933440, a_even=4802231769216, H_m=114569137164224
```

Pair-standard finite reductions include:

```text
n=5: lam_prime values 28,10,2; mu values 124,42,2
n=7: lam_prime values 7320,1460,84; mu values 13560,4100,132
```

## Decision

P15-S5 PASS as a reduction gate. The full odd theorem is not yet closed: `H_m` all-`m` certification and general odd pair-standard positivity remain open.

Next gate: P15-S6 `C_2 x C_2` square and closed forms.
