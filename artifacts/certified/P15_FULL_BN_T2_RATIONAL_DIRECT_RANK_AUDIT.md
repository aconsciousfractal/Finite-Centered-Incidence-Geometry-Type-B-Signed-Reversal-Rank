# P15 T2 Rational Direct-Rank Audit

Status: **PASS**

## Claim

```text
T_2^(8,3) direct Specht block ranks and R+ bounds have been recomputed exactly over Q; all selected listed defect CSV rows match the rational ranks.
```

## Rows

| lambda_minus | total dim | denom | rank_Q | R+ bound_Q | defect_Q | bucket_Q | listed CSV | CSV match |
|---|---:|---:|---:|---:|---:|---|---|---|
| `6` | 28 | 1 | 12 | 16 | 4 | EXTRA_DEFECT | True | True |
| `5.1` | 140 | 60 | 36 | 68 | 32 | EXTRA_DEFECT | True | True |
| `4.2` | 252 | 36 | 60 | 132 | 72 | EXTRA_DEFECT | True | True |
| `4.1.1` | 280 | 120 | 48 | 136 | 88 | EXTRA_DEFECT | True | True |
| `3.3` | 140 | 12 | 36 | 64 | 28 | EXTRA_DEFECT | True | True |
| `3.2.1` | 448 | 96 | 96 | 224 | 128 | EXTRA_DEFECT | True | True |
| `3.1.1.1` | 280 | 120 | 48 | 144 | 96 | EXTRA_DEFECT | True | True |
| `2.2.2` | 140 | 12 | 36 | 76 | 40 | EXTRA_DEFECT | True | True |
| `2.2.1.1` | 252 | 36 | 60 | 120 | 60 | EXTRA_DEFECT | True | True |
| `2.1.1.1.1` | 140 | 60 | 36 | 72 | 36 | EXTRA_DEFECT | True | True |
| `1.1.1.1.1.1` | 28 | 1 | 12 | 12 | 0 | EQUAL_R_BOUND | False | None |

## Method

The audit computes Young-seminormal Specht matrices over `Q`, takes a common denominator for the permutations that occur in the compressed fixed-plus operator and in the reversal operator, scales to integer matrices, and calls exact SymPy rank. For `R+`, the matrix ranked is `d I + R_scaled`, which has the same rank as `I + R` over `Q`.

## Sources

```text
scripts/p15_full_bn_t124_rational_direct_rank_audit.py
experiments/external_agents_2026_07_07/p15_full_bn_fingerprint_n8_kge2_defects.csv
```

## Boundary

Finite one-row fixed-plus rational direct-rank audit only; quotient-layer mechanism certificates remain separate, and no all-n full B_n fingerprint theorem is claimed.

Prossima task: Run the same rational direct-rank audit on the next family, in order T_4, T_2, then T_1; after all pass, update the P15 ledgers from modular to rational-direct where appropriate.
