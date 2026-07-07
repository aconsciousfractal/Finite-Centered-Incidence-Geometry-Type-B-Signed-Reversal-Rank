# P15 Full Bn Appendix Boundary Red-Team

Status: **BANKED / RED-TEAMED / APPENDIX-SCOPED**

Date: 2026-07-07

## Purpose

This note red-teams the finite full native `B_n` fingerprint material before it is used in a paper appendix.
It is not a new theorem gate for a full fingerprint classification. Its job is to decide what can be said safely, what remains open, and how the finite artifacts should be packaged for reviewers.

## Inputs Audited

```text
experiments/external_agents_2026_07_07/p15_full_bn_fingerprint_n8_kge2_defects.csv

docs/P15_FULL_BN_N8_ZERO_BLOCK_MECHANISMS_2026_07_07.md
docs/P15_FULL_BN_T1_K2_COMPLETE_FIXED_PLUS_ACCOUNTING_2026_07_07.md
docs/P15_FULL_BN_T2_K3_COMPLETE_FIXED_PLUS_ACCOUNTING_2026_07_07.md
docs/P15_FULL_BN_T4_K3_COMPLETE_FIXED_PLUS_ACCOUNTING_2026_07_07.md
docs/P15_FULL_BN_LAMBDA41_AND_T124_RATIONAL_UPGRADE_2026_07_07.md

certified/P15_FULL_BN_ZERO_FAMILY_GROUP_ALGEBRA_CERTIFICATE.md
certified/P15_FULL_BN_ZERO_FAMILY_CERTIFICATE.md
certified/P15_FULL_BN_PURE_PLUS_ZERO_CERTIFICATE.md
certified/P15_FULL_BN_PURE_PLUS_JOHNSON_DEGREE_DROP_CERTIFICATE.md
certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.md
certified/P15_FULL_BN_T1_RATIONAL_DIRECT_RANK_AUDIT.md
certified/P15_FULL_BN_T2_RATIONAL_DIRECT_RANK_AUDIT.md
certified/P15_FULL_BN_T4_RATIONAL_DIRECT_RANK_AUDIT.md
```

## Verdict

The appendix can safely say that the `n=8` full native `B_n` fingerprint scout refutes naive all-irrep saturation and that selected finite defect families are exactly accounted.
It must not say that P15 has a full native `B_n` fingerprint theorem or classification.

The two zero mechanisms are strong enough for finite appendix propositions:

1. `lambda_plus=(1,1)`, `lambda_minus=nu |- 6`: the block already vanishes in `Mat_28(Q[S_6])`, before applying any complement Specht module.
2. `lambda_plus=(2,2,2,1,1)`, `lambda_minus=empty`: the pure-plus zero has a finite conceptual explanation by Johnson degree-drop after sign twist, with an independent direct rational Specht zero check.

The extra-defect mechanisms are only partially accounted. They are useful as finite evidence and examples, but not as a classification.

## Exact Census After Current Closures

The external `n=8,k>=2` defect CSV has:

| slice | total fixed-plus families | total defect rows | accounted families | accounted rows | open families | open rows |
|---|---:|---:|---:|---:|---:|---:|
| `n=8,k=2` | 5 | 20 | 1 | 15 | 4 | 5 |
| `n=8,k=3` | 36 | 94 | 5 | 29 | 31 | 65 |
| total | 41 | 114 | 6 | 44 | 35 | 70 |

Accounted family-instances:

```text
n=8,k=2: lambda_plus=1
n=8,k=3: lambda_plus=1.1, 2, 4, 4.1, 2.2.2.1.1
```

Open family-instances include the remaining `n=8,k=2` families and the remaining `n=8,k=3` families, for example `lambda_plus=4.1.1` and beyond.

## Red-Team Findings

### Zero mechanism A: `(1.1 | nu)`, `nu |- 6`

Supported claim:

```text
For B~I_8(3,3), the fixed-plus block with lambda_plus=(1,1)
vanishes in Mat_28(Q[S_6]). Therefore all eleven Specht blocks
(lambda_plus=(1,1), lambda_minus=nu), nu |- 6, are zero.
```

Why this is strong:

```text
active_terms = 192
group_algebra_keys = 96
nonzero_group_algebra_keys = 0
cancelling_pairs = 96
bad_pairs = 0
signed_coefficient_patterns = {(-4, 4): 96}
```

This is cancellation before representation evaluation. It is not a modular accident and it is not tied to a specific `nu`.

Boundary:

```text
Finite n=8,k=3 proposition only.
No all-n lambda_plus=(1,1) vanishing theorem is claimed.
```

### Zero mechanism B: pure-plus `(2.2.2.1.1 | empty)`

Supported claim:

```text
M_(2.2.2.1.1 | empty)(B~I_8(3,3)) = 0.
```

The direct rational Specht certificate proves the block is exactly zero. The Johnson certificate gives the better finite explanation:

```text
(2,2,2,1,1)' = (5,3)
S^(2,2,2,1,1) ~= S^(5,3) tensor sgn
A_sgn(M^(3)) <= U_2
M^(3)/U_2 ~= S^(5,3)
```

So after sign twist, the operator lands in the pair-incidence submodule `U_2` of the 3-subset permutation module and kills the top quotient `S^(5,3)`.

Boundary:

```text
This is a finite Johnson/permutation-module proposition for n=8,k=3.
It is stronger than a black-box finite witness, but it is not an all-n Johnson theorem.
```

### Extra defects

Supported claim:

```text
The fixed-plus packages T_1^(8,2), T_2^(8,3), T_4^(8,3), and
lambda_plus=(4,1) at n=8,k=3 have finite accounting certificates.
The direct ranks and R+ bounds for T_1/T_2/T_4 have exact rational audits.
```

Boundary:

```text
The rational T_1/T_2/T_4 audit upgrades direct rank/bound evidence over Q.
It does not by itself convert every quotient-layer explanation into an all-n theorem.
The remaining 35 family-instances and 70 rows are still open.
```

## Paper Placement

Use this material only in an appendix or optional finite-atlas section.
The main theorem should remain the scoped two-standard-channel theorem, plus the separately scoped S12 separator theorem.
The full-fingerprint appendix should be framed as:

```text
A finite n=8 native B_n fingerprint boundary check and defect-atlas appendix.
```

Safe appendix propositions:

1. Exact finite cancellation for the eleven `(1.1|nu)` zero rows.
2. Exact finite Johnson degree-drop for the pure-plus zero row.
3. Exact finite accounting packages for the selected fixed-plus families already certified.
4. Explicit statement that the full fingerprint classification remains open.

Do not use:

```text
classification
full fingerprint theorem
all irreducible saturation theorem
all-n Johnson theorem
native B_n theorem beyond the two standard channels
```

## Public-Repo Consequence

The public reviewer package should separate replay tiers:

1. Core theorem replay and manuscript certificates.
2. S12 separator replay.
3. Optional full-fingerprint appendix replay, with dependencies declared explicitly.

The full-fingerprint scripts are useful, but they are heavier than the Type-A public replay style because several use `numpy` and the rational direct-rank audit uses `sympy`.
Therefore the public repo should either declare these dependencies in `requirements.txt` or put the full-fingerprint replay in an optional appendix tier.

## Public Status

The finite full-fingerprint boundary appendix is integrated as an appendix-scoped section, and the public repository follows the Type-A-style replay/package structure.
