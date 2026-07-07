# P15-S11A-H True-Scalar Edge-Insertion Verifier

Schema: `p15.s11a.h_true_scalar_edge_insertion.v1`
Gate: P15-S11A-H true-scalar edge-insertion verifier
Status: **PASS**

## Decision

S11A-H passes: true center-nonfixed edge-insertion local algebra is certified, small edge rows are exact, and the all-m Phi-kernel is proven by moment integration by parts.

## Checks

- `free_pair_constant_transfer_Az`: **PASS**
  Coefficient [u^0 v^0] of a full free reversal pair is Az=1+(2z-6)s+(z^2-4z+5)s^2.
- `forced_edge_class_formulas`: **PASS**
  The center split gives one fixed target, one reversal target, and 2m-2 ordinary targets per source; local factors are singleton fixed, singleton reversal, and two hooks.
- `direct_small_edge_rows_m2_m4`: **PASS**
  Direct exact edge rows close the m=2,3,4 cases outside the A^(m-5) kernel range.
- `edge_insertion_vs_agent1_phi_smoke`: **PASS**
  Diagnostic exact rows for the complete local-vs-compact equality; proof is supplied by the small rows plus moment reduction.
- `compact_difference_kernel_residual_Q`: **PASS**
  Over Z[m,s], E_z(local-(2m-2)R_m)=A(s)^(m-5)Q_m(s) for m>=5.
- `phi_kernel_moment_reduction`: **PASS**
  Integration-by-parts moment recurrence reduces Phi_{2m-1}(A^(m-5)Q_m) to zero for all m>=5.

## Kernel Proof

The exact residual is

```text
E_z(local_edge_insertion - (2m-2)R_m) = A(s)^(m-5) Q_m(s),
E_z(F)=F(1,s)-d_zF(1,s).
```

For `m>=5`, the substitution `s=1/(2x)` reduces the kernel to moments `I_j=int_0^infty exp(-x)B(x)^(m-5)x^j dx`, `B=2x^2-4x+1`. Integration by parts gives

```text
j I_{j-1}+(-4j-4m+15)I_j+(2j+4m-12)I_{j+1}-2I_{j+2}=-delta_{j0}.
```

The verifier reduces the transformed target to zero in the basis `I_0`, `I_1`, `1`. The small rows `m=2,3,4` are checked directly.

## Q_m Terms

```text
m^0 s^0: 2
m^1 s^0: -2
m^0 s^1: -26
m^1 s^1: 22
m^2 s^1: 4
m^0 s^2: 108
m^1 s^2: 8
m^2 s^2: -148
m^3 s^2: 32
m^0 s^3: -420
m^1 s^3: 32
m^2 s^3: 420
m^3 s^3: 32
m^4 s^3: -64
m^0 s^4: 832
m^1 s^4: -212
m^2 s^4: -172
m^3 s^4: -896
m^4 s^4: 448
m^0 s^5: -904
m^1 s^5: 432
m^2 s^5: -1128
m^3 s^5: 2816
m^4 s^5: -1216
m^0 s^6: 816
m^1 s^6: -1072
m^2 s^6: 2624
m^3 s^6: -3968
m^4 s^6: 1600
m^0 s^7: -816
m^1 s^7: 1536
m^2 s^7: -2320
m^3 s^7: 2624
m^4 s^7: -1024
m^0 s^8: 416
m^1 s^8: -560
m^2 s^8: 400
m^3 s^8: -512
m^4 s^8: 256
m^0 s^9: 64
m^1 s^9: -288
m^2 s^9: 352
m^3 s^9: -128
```
