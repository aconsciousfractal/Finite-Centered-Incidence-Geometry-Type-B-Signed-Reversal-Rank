# P15-S11A Scalar-Link Symbolic Transfer Verifier

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S11A scalar-link symbolic transfer verifier
Date: 2026-07-07
Status: PASS for transfer algebra only

## Decision

The symbolic transfer identities pass over `Z[m,t]`. This verifies that the S11A intermediate coefficient-extraction formulas reduce algebraically to the displayed univariate targets `P_m(t)` and `L_m(t)`.

This local transfer certificate does not by itself close S11A, but downstream S11A-H and S11A-K certificates now close both all-`m` true-scalar links.

## Checks

```text
H_transfer_agent1_to_univariate_P_m  PASS  residual terms: 0
K_bivariate_to_univariate_L_m        PASS  residual terms: 0
```

For `H_m`, the verifier checks the local identity obtained from

```text
R_m(z,t)=az(z,t)^(m-2) (1+(z-2)t)(1-t)(1-2t),
S_m=Phi(R_m(1,t)-d_z R_m(1,t)).
```

After factoring `A(t)^(m-3)`, the local polynomial equals the displayed `(1-t)P_m(t)` target.

For `K_m`, the verifier starts from the S9 bivariate target

```text
c * ( a^(m-1)*((1+u)sx+(1+v)sy) - 4*a^(m-2)*h^2 )
```

extracts `[u v]`, factors `B(t)^(m-4)`, and verifies that the remaining local polynomial is exactly the displayed degree-eight `L_m(t)`.

## Boundary

No held-out range was extended. No symbolic search was performed. The verifier checks identities only. The former downstream closure tasks are now closed:

```text
S11A-H: closed by `P15_S11A_H_TRUE_SCALAR_EDGE_INSERTION_VERIFIER` plus S8.
S11A-K: closed by `P15_S11A_K_TRUE_SCALAR_BIVARIATE_VERIFIER`, `P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE`, and S9C.
```

## Artifacts

- `scripts/p15_s11a_scalar_link_symbolic_verifier.py`
- `certified/P15_S11A_SCALAR_LINK_SYMBOLIC_VERIFIER.json`
