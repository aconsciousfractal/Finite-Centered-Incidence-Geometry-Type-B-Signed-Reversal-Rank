# P15 Full Bn n=8 Zero-Block Mechanisms

Status: **BANKED / EXACT FINITE CERTIFICATES PASS / PURE-PLUS JOHNSON MECHANISM VERIFIED / APPENDIX RED-TEAMED**

Date: 2026-07-07

## Purpose

This note isolates the exact zero rows in the full native `B_n` fingerprint scout for `B~I_8(3,3)`.
It is a finite defect-atlas artifact. It does not promote a full native `B_n` fingerprint theorem.

## Zero Census

The n=8, k=3 atlas has exactly 12 `ZERO` rows.
They split into two mechanisms:

| mechanism | zero rows | certificate | reading |
|---|---:|---|---|
| `(1.1 | nu)`, `nu |- 6` | 11 | `P15_FULL_BN_ZERO_FAMILY_GROUP_ALGEBRA_CERTIFICATE` | uniform pair-cancellation in `Q[S_6]` before applying Specht modules |
| `(2.2.2.1.1 | empty)` | 1 | `P15_FULL_BN_PURE_PLUS_JOHNSON_DEGREE_DROP_CERTIFICATE`; independently `P15_FULL_BN_PURE_PLUS_ZERO_CERTIFICATE` | Johnson degree-drop after sign twist; independently checked by rational `S_8` Specht sum |

So the original 12 zero rows from the modular atlas are now all accounted for by exact local certificates.

## Mechanism A: `(1.1 | nu)`, `nu |- 6`

For `lambda_plus=(1,1)`, the plus-sector Specht module is the sign representation of `S_2`.
The layer operator for `B~I_8(3,3)` can be computed before applying a complement Specht module:

```text
T_{1.1} in Mat_{28}(Q[S_6]).
```

The certificate groups active signed terms by source plus subset, destination plus subset, and complement permutation in `S_6`.
The exact audit is:

```text
active_terms = 192
group_algebra_keys = 96
nonzero_group_algebra_keys = 0
cancelling_pairs = 96
bad_pairs = 0
signed_coefficient_patterns = {(-4, 4): 96}
```

For every fixed complement key, the two surviving completions have the same complement map, swap the two plus images, have the same unsigned coefficient `4`, and have opposite `S_2` sign. Thus each key contributes `4 sigma - 4 sigma = 0` in `Q[S_6]`. Applying any Specht module `S^nu`, `nu |- 6`, preserves zero.

A finite scout through diagonal rows `4 <= n <= 8` found this total pair-cancellation only at `(n,k)=(8,3)`. Therefore this is not being promoted as a generic all-`n` `lambda_plus=(1,1)` vanishing theorem.

## Mechanism B: `(2.2.2.1.1 | empty)`

The remaining zero row is not part of the `(1.1|nu)` family. It is the pure-plus block

```text
M_{(2.2.2.1.1 | empty)}(B~I_8(3,3)) = 0.
```

The direct seminormal certificate is exact, but the better explanation is a Johnson degree-drop.
Since

```text
(2,2,2,1,1)' = (5,3),
S^(2,2,2,1,1) ~= S^(5,3) tensor sgn,
```

the pure-plus block can be sign-twisted and tested on the 3-subset permutation module `M^(3)`.
Let `A_sgn=sum_pi w(pi) sign(pi) pi`, where `w(pi)` is the number of sign choices producing three positive identity hits and three positive reversal hits.
The Johnson certificate proves

```text
A_sgn(M^(3)) <= U_2,
```

where `U_2` is the pair-incidence submodule spanned by the vectors of 3-subsets containing a fixed pair. Therefore `A_sgn` kills

```text
M^(3)/U_2 ~= S^(5,3),
```

and untwisting by sign gives the original pure-plus zero.

Exact Johnson audit:

```text
nonzero_underlying_permutations = 102
weight_patterns = {(3,3,4): 96, (4,4,16): 6}
rank_pair_incidence_U2 = 28
rank_sign_twisted_operator = 16
formula_residual_nonzero_entries = 0
```

The verified coefficient formula is:

```text
A_sgn e_T = sum_{|P|=2} c(P,T) v_P,

c({x,y},T) = -8 * (
    1[x in T and y in rho(T)]
  + 1[y in T and x in rho(T)]
  - 1[{x,y} is a rho-pair] * 1[|{x,y} cap T| = 1]
).
```

The previous rational `S_8` Specht sum remains an independent exact check:

```text
specht_dim = 28
common_denominator_seen = 2160
nonzero_entries = 0
```

This is a finite `n=8,k=3` mechanism, not an all-`n` theorem.

## Analogy Notes

P12 already uses a quotient of a 2-subset permutation module plus sign twists as an independent rational verifier. P05/P04 contain Johnson/fixed-weight language but keep theorem-level Johnson claims source-locked. Therefore the safe P15 phrasing is a local finite Johnson/permutation-module certificate, not a general Johnson-scheme theorem.

## Sources

```text
scripts/p15_full_bn_zero_block_certificate.py
scripts/p15_full_bn_zero_family_certificate.py
scripts/p15_full_bn_zero_family_group_algebra_certificate.py
scripts/p15_full_bn_lambda11_group_algebra_scout.py
scripts/p15_full_bn_pure_plus_zero_certificate.py
scripts/p15_full_bn_pure_plus_johnson_degree_drop.py

artifacts/certified/P15_FULL_BN_ZERO_BLOCK_CERTIFICATE.md
artifacts/certified/P15_FULL_BN_ZERO_FAMILY_CERTIFICATE.md
artifacts/certified/P15_FULL_BN_ZERO_FAMILY_GROUP_ALGEBRA_CERTIFICATE.md
artifacts/certified/P15_LAMBDA11_GROUP_ALGEBRA_SCOUT_N8.md
artifacts/certified/P15_FULL_BN_PURE_PLUS_ZERO_CERTIFICATE.md
artifacts/certified/P15_FULL_BN_PURE_PLUS_JOHNSON_DEGREE_DROP_CERTIFICATE.md
```

## Red-Team 2026-07-07

The two zero mechanisms survive appendix-boundary red-team.
The `(1.1|nu)` family is a finite pre-Specht group-algebra cancellation in `Mat_28(Q[S_6])`, so it is stronger than eleven separate modular zero checks.
The pure-plus row is not merely a black-box finite witness: after sign twist it is a finite Johnson degree-drop proposition, independently checked by a direct rational Specht sum.

The safe paper use is still finite and local: these are `n=8,k=3` mechanisms inside the defect atlas, not a full native `B_n` fingerprint theorem.

## Boundary

This closes the exact finite zero-row accounting for the `n=8,k=3` atlas and gives conceptual finite mechanisms for both zero families.
It does not classify the remaining `EXTRA_DEFECT` rows and does not prove a full native `B_n` fingerprint theorem.

Current defect-atlas accounting after the latest fixed-plus closures:

```text
n=8,k=2: 5 defect families, 20 rows; accounted 1 family / 15 rows; open 4 families / 5 rows.
n=8,k=3: 36 defect families, 94 rows; accounted 5 families / 29 rows; open 31 families / 65 rows.
total:    41 defect families, 114 rows; accounted 6 families / 44 rows; open 35 families / 70 rows.
```

The `n=8,k=2` `lambda_plus=1` family is now accounted elsewhere, but the other `n=8,k=2` families remain open.
The rational `T_1/T_2/T_4` audits upgrade direct rank/bound evidence over `Q`; they do not turn the finite atlas into an all-`n` classification.

Public status: the finite full-fingerprint appendix text and public reproduce package are integrated, with this material kept separate from the main theorem.
