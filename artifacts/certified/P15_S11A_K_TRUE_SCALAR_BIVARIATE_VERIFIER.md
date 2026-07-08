# P15-S11A-K True-Scalar To Bivariate Verifier

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S11A-K true-scalar to bivariate local verifier
Date: 2026-07-07
Status: PASS for local signed-rook factors; global product/Phi closed downstream

## Decision

The local signed-rook factors for the true odd pair-standard `k=1` scalar match the S9 bivariate target. Together with `P15_S11A_SCALAR_LINK_SYMBOLIC_VERIFIER`, this gives a symbolic proof route for S11A-K:

```text
true lambda_prime entry formula
  -> S9 bivariate target
  -> displayed L_m(t) univariate target
  -> S9C recurrence/positivity certificate.
```

This certificate does not rely on held-out values and does not search for formulas.

## True Scalar Formula

For odd `n=2m+1`, choose a non-center reversal pair `{i,rho(i)}` and a generic value `j` outside that pair and the center. S5 records

```text
lambda_prime = M_pair[i,i] + M_pair[rho(i),i] - 2*M_pair[j,i].
```

The verifier checks the local signed-rook factors that convert this true entry formula into the S9 bivariate expression.

## Checks

```text
free_reversal_pair_transfer_a              PASS  residual terms: 0
odd_center_transfer_c                      PASS  residual terms: 0
same_pair_forced_Mii_plus_Mrhoi_i          PASS  residual terms: 0
generic_forced_minus_2_Mji_term            PASS  residual terms: 0
lambda_prime_true_entry_shape              PASS  residual terms: 0
```

The global formula certified by the local factors is

```text
lambda_prime(2m+1,1)
 = [u v] Phi_{2m}( c*(a^(m-1)*same - a^(m-2)*generic) ),
generic = 4*h^2.
```

## Boundary

This is a local algebra certificate. The standard product over independent non-center reversal pairs and the `Phi_{2m}` factorial functional are now closed by `P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE`. S11A-H is closed separately.

## Artifacts

- `scripts/p15_s11a_k_true_scalar_bivariate_verifier.py`
- `artifacts/certified/P15_S11A_K_TRUE_SCALAR_BIVARIATE_VERIFIER.json`
