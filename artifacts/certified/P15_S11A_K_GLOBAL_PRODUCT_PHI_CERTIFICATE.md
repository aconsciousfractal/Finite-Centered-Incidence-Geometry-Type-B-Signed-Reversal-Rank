# P15-S11A-K Global Product/Phi Certificate

Schema: `p15.s11a.k_global_product_phi.v1`
Gate: P15-S11A-K global product/Phi certificate
Status: **PASS**

## Decision

S11A-K global product/Phi assembly is certified: the true lambda_prime entry formula gives the S9 bivariate target for all m>=2; together with transfer and S9C this closes the pair-standard k=1 scalar link.

## Formula

For odd `n=2m+1`, `m>=2`, choose a non-center reversal pair `{i,rho(i)}`, the center `c0`, and a generic source `j` outside that pair and the center. S5 gives

```text
lambda_prime = M_pair[i,i] + M_pair[rho(i),i] - 2*M_pair[j,i].
```

With `X=u-1`, `Y=v-1`, `Z=uv-1`, define

```text
a = 1+2t(X+Y)+t^2(X^2+Y^2)
c = 1+tZ
sx = 1+tX
sy = 1+tY
h = 1+t(X+Y)
same = (1+u)sx+(1+v)sy
generic = 4h^2
```

Then the all-`m` product/Phi identity is

```text
lambda_prime(2m+1,1) = [uv] Phi_{2m}( c*(a^(m-1)*same - a^(m-2)*generic) ).
```

## Proof

A single matrix entry in `lambda_prime` fixes one source-target assignment in an odd board of size `2m+1`. The remaining completion board has `2m` rows and `2m` columns.

For `M_pair[i,i]+M_pair[rho(i),i]`, the forced edge lies in the distinguished reversal pair. The center contributes `c`, the unused mate in the distinguished pair contributes `sx` or `sy`, and the other `m-1` non-center reversal pairs contribute the free factor `a`. The two forced signs have weights `1+u` and `1+v`, so the same-pair term is `c*a^(m-1)*same`.

For `M_pair[j,i]`, the forced ordinary edge touches two non-center reversal pairs. Those pairs leave two independent hook factors, giving `h^2`. The forced ordinary sign weight is `2`, and the coefficient `-2` in `lambda_prime` gives the subtracted magnitude `4h^2=generic`. The center contributes `c`, and the remaining `m-2` pairs contribute `a^(m-2)`. Thus the generic term is subtracted as `c*a^(m-2)*generic`.

The factors use disjoint row and column sets, so the signed-rook polynomial of the union is the product of the block polynomials. If a term selects `r` nonattacking special cells in the remaining `2m` by `2m` board, the unmatched rows and columns can be completed by an arbitrary signed bijection, giving `2^(2m-r)(2m-r)!` completions. This is exactly `Phi_{2m}`.

The case `m=1` has no pair-standard `lambda_prime` eigenspace. For `m=2,3`, the bivariate product has exponents `m-1` and `m-2` nonnegative, while the univariate `B(t)^(m-4)L_m(t)` form is not invoked; S9C handles these rows directly. For `m>=4`, the S11A transfer verifier proves the bivariate target collapses to `B(t)^(m-4)L_m(t)`.

## Checks

- `dependency_verifiers_pass`: **PASS**
- `same_pair_product_accounting`: **PASS**
- `generic_pair_product_accounting`: **PASS**
- `phi_2m_completion_functional`: **PASS**
- `edge_policy_m1_m2_m3`: **PASS**

## Chain

- P15_S11A_K_TRUE_SCALAR_BIVARIATE_VERIFIER: local true entry factors pass.
- P15_S11A_SCALAR_LINK_SYMBOLIC_VERIFIER: S9 bivariate target collapses to B(t)^(m-4)L_m(t).
- P15_S9C_KM_RECURRENCE_CERTIFICATE: displayed K_m formula sequence is positive for all m>=2.

## Boundary

This certificate closes the S11A-K formula-to-true-scalar link. It does not broaden the public theorem beyond the S10 two-standard-channel boundary, and it makes no Type-D, full-fingerprint, separator, or enumeration-novelty claim.
