# P15 lambda_plus=4.1 k=3 Plus-Standard Certificate

Status: **PASS**

## Claim

```text
The lambda_plus=(4,1), n=8,k=3 fixed-plus family is finitely accounted: all three complement rows are the plus-standard quotient of the 5-point plus permutation module modulo constants, and their ranks match the defect CSV.
```

## Profiles

| lambda_minus | natural rank | lower rank | quotient rank | direct Specht rank | R+ bound | defect | bucket | CSV match |
|---|---:|---:|---:|---:|---:|---:|---|---|
| `3` | 104 | 56 | 76 | 76 | 112 | 36 | EXTRA_DEFECT | True |
| `2.1` | 232 | 112 | 176 | 176 | 224 | 48 | EXTRA_DEFECT | True |
| `1.1.1` | 128 | 56 | 100 | 100 | 112 | 12 | EXTRA_DEFECT | True |

## Mechanism

For `lambda_plus=(4,1)`, the plus-side Specht module is the standard quotient of the 5-point permutation module by constants. The certificate builds that natural module over each plus subset fiber, checks that the layer operator maps fiber constants back into fiber constants, and then compares the quotient image rank with the direct `S^(4,1) tensor S^nu` Specht block rank.

## Checks

### `3`

```json
{
  "lower_rank_is_expected_fiber_constants": true,
  "source_lower_maps_inside_lower": true,
  "quotient_rank_matches_direct_specht_rank": true,
  "direct_minus_dim_matches_natural_minus_dim": true,
  "csv_match": true,
  "extra_defect": true
}
```

### `2.1`

```json
{
  "lower_rank_is_expected_fiber_constants": true,
  "source_lower_maps_inside_lower": true,
  "quotient_rank_matches_direct_specht_rank": true,
  "direct_minus_dim_matches_natural_minus_dim": true,
  "csv_match": true,
  "extra_defect": true
}
```

### `1.1.1`

```json
{
  "lower_rank_is_expected_fiber_constants": true,
  "source_lower_maps_inside_lower": true,
  "quotient_rank_matches_direct_specht_rank": true,
  "direct_minus_dim_matches_natural_minus_dim": true,
  "csv_match": true,
  "extra_defect": true
}
```

## Sources

```text
scripts/p15_full_bn_lambda41_k3_plus_standard_certificate.py
artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.json
artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.md
experiments/external_agents_2026_07_07/p15_full_bn_fingerprint_n8_kge2_defects.csv
```

## Boundary

Finite n=8,k=3 lambda_plus=(4,1) mechanism certificate only; no all-n fixed-plus theorem and no full native B_n fingerprint theorem.

Prossima task: Optional rational direct-rank audits for T_1/T_2/T_4 are separate replay-tier checks; no full native B_n fingerprint classification is promoted.
