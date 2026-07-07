# P15 T1 Rational Direct-Rank Audit

Status: **PASS**

## Claim

```text
T_1^(8,2) direct Specht block ranks and R+ bounds have been recomputed exactly over Q; all selected listed defect CSV rows match the rational ranks.
```

## Rows

| lambda_minus | total dim | denom | rank_Q | R+ bound_Q | defect_Q | bucket_Q | listed CSV | CSV match |
|---|---:|---:|---:|---:|---:|---|---|---|
| `7` | 8 | 1 | 0 | 4 | 4 | ZERO | True | True |
| `6.1` | 48 | 60 | 12 | 24 | 12 | EXTRA_DEFECT | True | True |
| `5.2` | 112 | 720 | 36 | 56 | 20 | EXTRA_DEFECT | True | True |
| `5.1.1` | 120 | 360 | 48 | 60 | 12 | EXTRA_DEFECT | True | True |
| `4.3` | 112 | 72 | 36 | 56 | 20 | EXTRA_DEFECT | True | True |
| `4.2.1` | 280 | 1440 | 108 | 140 | 32 | EXTRA_DEFECT | True | True |
| `4.1.1.1` | 160 | 720 | 72 | 80 | 8 | EXTRA_DEFECT | True | True |
| `3.3.1` | 168 | 288 | 60 | 84 | 24 | EXTRA_DEFECT | True | True |
| `3.2.2` | 168 | 288 | 60 | 84 | 24 | EXTRA_DEFECT | True | True |
| `3.2.1.1` | 280 | 1440 | 108 | 140 | 32 | EXTRA_DEFECT | True | True |
| `3.1.1.1.1` | 120 | 360 | 48 | 60 | 12 | EXTRA_DEFECT | True | True |
| `2.2.2.1` | 112 | 72 | 36 | 56 | 20 | EXTRA_DEFECT | True | True |
| `2.2.1.1.1` | 112 | 720 | 36 | 56 | 20 | EXTRA_DEFECT | True | True |
| `2.1.1.1.1.1` | 48 | 60 | 12 | 24 | 12 | EXTRA_DEFECT | True | True |
| `1.1.1.1.1.1.1` | 8 | 1 | 0 | 4 | 4 | ZERO | True | True |

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
