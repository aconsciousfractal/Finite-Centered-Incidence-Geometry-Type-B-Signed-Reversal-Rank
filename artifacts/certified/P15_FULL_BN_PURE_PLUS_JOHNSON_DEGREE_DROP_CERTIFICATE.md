# P15 Pure-Plus Johnson Degree-Drop Certificate

Status: **PASS**

## Claim

```text
The pure-plus zero M_{(2.2.2.1.1|empty)}(B~I_8(3,3)) is explained by Johnson degree-drop after sign twist: A_sgn(M^(3)) is contained in the pair-incidence submodule U_2, so the top quotient S^(5,3) is killed.
```

## Mechanism

The pure-plus zero block is `lambda=(2,2,2,1,1)`. Its conjugate is `(5,3)`, so

```text
S^(2,2,2,1,1) ~= S^(5,3) tensor sgn.
```

After sign twist, the target is the top quotient `S^(5,3)` of the permutation module on 3-subsets of `[8]`.
Let `v_P` be the incidence vector of all 3-subsets containing the pair `P`. The pair-incidence span `U_2=<v_P>` has rank 28 and contains only the lower Johnson layers.
The certificate verifies exactly that the sign-twisted weighted operator maps every 3-subset basis vector into `U_2`.

## Coefficient Formula

For a pair `P={x,y}` and a 3-subset `T`, with `rho(i)=7-i`, the verified formula is

```text
c({x,y},T) = -8 * (
    1[x in T and y in rho(T)]
  + 1[y in T and x in rho(T)]
  - 1[{x,y} is a rho-pair] * 1[|{x,y} cap T| = 1]
)
A_sgn e_T = sum_{|P|=2} c(P,T) v_P.
```

Thus `A_sgn(M^(3)) <= U_2`, so `A_sgn` kills `M^(3)/U_2 ~= S^(5,3)`. Untwisting by sign gives the pure-plus zero `S^(2,2,2,1,1)`.

## Exact Audit

| field | value |
|---|---:|
| `subset_module_dim` | 56 |
| `pair_count` | 28 |
| `nonzero_underlying_permutations` | 102 |
| `layer_size` | 480 |
| `rank_pair_incidence_U2` | 28 |
| `rank_sign_twisted_operator` | 16 |
| `formula_residual_nonzero_entries` | 0 |
| `formula_residual_max_abs` | 0 |

## Weight Patterns

```text
{'(3, 3, 4)': 96, '(4, 4, 16)': 6}
```

## Boundary

Finite n=8,k=3 pure-plus mechanism only; no all-n full native B_n fingerprint theorem or classification theorem.

Prossima task: test whether the two `n=8,k=2` zero rows also admit a pre-Specht or Johnson/permutation-module degree-drop explanation.
