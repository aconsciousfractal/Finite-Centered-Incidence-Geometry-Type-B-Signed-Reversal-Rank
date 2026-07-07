# P15-S11 Manuscript Lemmas

Date: 2026-07-07
Status: INTEGRATED into the split public manuscript. These lemmas summarize the S4/S5/S8/S9/S9C/S11 certificate chain; they do not broaden the theorem.

## Red-Team Status

S11B positivity, S11C nonempty rows, and S11D raw-rank smoke checks have survived the current red-team round. S11A-H is now closed: the true odd `V_ref`, `k=1` scalar link is certified by local edge-insertion algebra, direct small edge rows, and an all-`m` Phi-kernel proof by moment integration by parts.

The S11A-K global product/Phi certificate now supplies the missing standard product and `Phi_{2m}` argument, and post-fix red-team found no blocking issue. Thus both S11A-H and S11A-K are closed internally; S11B/S11C remain manuscript-polish items, not proof blockers.

## Lemma S11A-H: `H_m` formula target for the odd `V_ref`, `k=1`, scalar residual

For `m>=2`, write `n=2m+1`. The intended odd reference-channel scalar link is

```text
a(2m+1,1) = a(2m,1) + H_m.
```

The edge value is `H_2=22`. For `m>=3`, put

```text
A(t)=1-4t+2t^2,
Phi_N(sum_r h_r t^r)=sum_r h_r 2^(N-r)(N-r)!,
P_m(t)=(2m-2)(2t-1)(2m t^3-(4m+2)t^2+(2m+2)t-1)-A(t)^2.
```

Then the formula target is

```text
H_m = Phi_{2m-1}((1-t) A(t)^(m-3) P_m(t)).
```

Proof outline. S5 gives the odd split `M_ref=2a P_+^{pairs}+a_c E_center`. In the `k=1` row the trace relation reduces the non-center scalar to the even scalar plus the center-nonfixed correction:

```text
(2m) a(2m+1,1) = (2m) a(2m,1) + C_m.
```

Thus it is enough to identify `C_m=2m H_m`. The center-nonfixed terms are intended to be counted by inserting the odd center into one functional edge of an even `2m` non-center configuration. The two local insertion alternatives are encoded by the Ryser/factorial functional with base `A(t)`. Their difference is

```text
H_m=(2m-2)L_m-B_m,
B_m=Phi_{2m-1}((1-t)A(t)^(m-1)),
L_m=Phi_{2m-1}((1-t)(2t-1)A(t)^(m-3)
      (2m t^3-(4m+2)t^2+(2m+2)t-1)),
```

which combines to the stated polynomial `P_m(t)`. The S11A-H edge-insertion verifier derives the true center-nonfixed scalar contribution from forced even-board edges (fixed, reversal, ordinary) and checks the local signed-rook factors. It reduces the comparison with the compact Agent-1 `R_m` expression to the explicit residual `Phi_{2m-1}(A(s)^(m-5)Q_m(s))`; direct exact rows close `m=2,3,4`, and the moment integration-by-parts recurrence proves the residual is zero for every `m>=5`. The S11A symbolic transfer verifier checks the subsequent local reduction to the displayed `P_m(t)` target over `Z[m,t]`, and S8 proves positivity of the formula sequence `H_m` for all `m>=2`. Thus S11A-H is proof-complete internally; finite smoke rows are diagnostics only.

## Lemma S11A-K: `K_m` formula target for the odd pair-standard `k=1` scalar

For odd `n=2m+1`, the pair-standard channel has the eigenvalue `lambda_prime` on the non-center reversal-plus standard subspace. For `m>=2`, set

```text
K_m=(m-1)lambda_prime(2m+1,1).
```

The edge values are

```text
K_2=28,  K_3=14640.
```

For `m>=4`, define

```text
B(t)=2t^2-4t+1,
Phi_N(sum_r h_r t^r)=sum_r h_r 2^(N-r)(N-r)!,
```

and

```text
L_m(t)=2(2m-3)t
     +2(-4m^2-8m+19)t^2
     +2(24m^2+18m-68)t^3
     +2(-56m^2-24m+122)t^4
     +2(56m^2+24m-92)t^5
     +2(-12m^2-44m+12)t^6
     +2(-16m^2+56m)t^7
     +2(8m^2-24m+8)t^8.
```

Then the formula target is

```text
lambda_prime(2m+1,1)=Phi_{2m}(B(t)^(m-4)L_m(t)),
K_m=(m-1)Phi_{2m}(B(t)^(m-4)L_m(t)).
```

The S11A-K global product/Phi certificate checks the remaining standard product over independent pairs and the `Phi_{2m}` functional: same-pair blocks give `c*a^(m-1)*same`, generic blocks give the subtracted `c*a^(m-2)*generic` with `generic=4h^2`, and a forced entry leaves `2m` rows and columns, so `Phi_{2m}` counts completions. Therefore S11A-K is proof-complete internally.

## Lemma S11B: positivity from the certified order-three recurrences

Let `U_m` be one of the formula sequences `H_m` or `K_m`, and suppose it satisfies a certified recurrence

```text
p_0(m)U_m+p_1(m)U_{m+1}+p_2(m)U_{m+2}+p_3(m)U_{m+3}=0.
```

Set `rho_m=U_{m+1}/U_m`. Whenever the base rows have `U_m>0` and the needed adjacent ratios `rho_m>=4m^2`, the recurrence gives

```text
rho_{m+2}=(-p_2(m)-p_1(m)/rho_{m+1}-p_0(m)/(rho_m rho_{m+1}))/p_3(m).
```

For `H_m`, S8 verifies for all `m>=2` that `p_0>0`, `p_1>0`, `-p_2>0`, `p_3>0`, and

```text
D_H(m)=16m^2(m+1)^2(-p_2)-4m^2p_1-p_0
       -64m^2(m+1)^2(m+2)^2p_3 > 0.
```

Together with the base rows `m=2,3,4`, where S8 checks both `H_m>0` and `rho_m>=4m^2`, this proves `rho_m>=4m^2` and `H_m>0` for every `m>=2`.

For `K_m`, S9C verifies the base rows `m=2..8`, including `K_m>0` and `rho_m>=4m^2` for those `m`, and for `m>=7` verifies `p_0>0`, `p_1<=0`, `-p_2>0`, `p_3>0`, and

```text
D_K(m)=16m^2(m+1)^2(-p_2)-p_0
       -64m^2(m+1)^2(m+2)^2p_3 > 0.
```

Here the nonpositive `p_1` term only improves the lower bound. Therefore the adjacent base ratios `rho_7` and `rho_8` seed the recurrence from `m=7`, `rho_m>=4m^2` propagates onward, and `K_m>0` for every `m>=2`.

## Lemma S11C: nonempty positive diagonal rows

For the theorem range `k>=1`, the diagonal layer `B~I_n(k,k)` is nonempty exactly in the following cases.

For `n=2m`:

```text
1 <= k < m, or k=m with m even.
```

For `n=2m+1`:

```text
1 <= k <= m, or k=m+1 with m even.
```

The row `k=0` is nonempty for every `n>=2`, but it is outside the P15 theorem statement.

Proof idea. Prescribe two sets of positive-hit positions: `F` for fixed hits and `R` for reversal hits. Off the odd center, a valid partial prescription must have no domain receiving two different images and must have distinct prescribed images. In particular, a non-center position cannot simultaneously be a positive fixed hit and a positive reversal hit.

For `n=2m`, the upper bound is `k<=m` because the `k` fixed-hit domains and `k` reversal-hit domains are disjoint. For `k<m`, construct examples by reversal-pair blocks. If `k=2a`, choose `a` reversal pairs to be fully fixed-positive and `a` disjoint reversal pairs to be fully reversal-positive. If `k=2a+1`, do the same for `2a` hits and use two further reversal pairs: one supplies a single positive fixed hit and the other supplies a single positive reversal hit. This uses `k` pairs when `k` is even and `k+1` pairs when `k` is odd, so it is possible for every `k<m`. At the top row `k=m`, all non-center positions are prescribed. If `I` is the set of positions sent as fixed, the complement is sent as reversed, and injectivity forces `I=rho(I)`. Thus `I` is a union of reversal pairs and has even cardinality; the top row exists exactly when `m` is even.

For `n=2m+1`, the center is fixed and reversed at the same time. Give it positive sign to contribute one fixed hit and one reversal hit. If `k=1`, use the empty non-center prescription; if `k>1`, apply the even `2m` construction to the non-center positions with `k-1` hits of each kind. This gives exactly `1<=k<=m`, plus the top row `k=m+1` when `m` is even.

The bounded S11D smoke check verifies this positive-row predicate through `n=6`, and S5 records the same odd top-row fence in its center-scalar sweep.

## Integration Decision

S11B, S11C, S11D, S11A-H, and S11A-K are red-team coherent. S11A-K is closed by local true-scalar factors, global product/Phi, bivariate-to-`L_m(t)` transfer, and S9C positivity. The theorem statement remains the S10-scoped statement only: nonempty diagonal `B~I_n(k,k)`, two ambient standard channels, no Type-D/general Weyl/full-fingerprint/enumeration-novelty claim.
