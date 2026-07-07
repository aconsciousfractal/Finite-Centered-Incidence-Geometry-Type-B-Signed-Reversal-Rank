# P15 Full Bn Zero Family Group-Algebra Certificate

Status: **PASS**

## Claim

```text
The lambda_plus=(1,1) block of B~I_8(3,3) is zero in every complement Q[S_6] block before applying Specht modules.
```

## Interpretation

For `lambda_plus=(1,1)`, the plus-sector Specht module is the sign representation of `S_2`.
After summing signs over `B~I_8(3,3)`, every surviving complement `S_6` basis term occurs exactly twice: the two terms have the same complement map, the two plus images swapped, equal unsigned coefficient `4`, and opposite plus-sector sign.
Therefore the alternating sum is already zero in the relevant block of `Q[S_6]`. Applying any Specht module `S^nu`, `nu |- 6`, preserves zero.

## Audit

| field | value |
|---|---:|
| `subset_count` | 28 |
| `permutations_scanned` | 40320 |
| `active_terms` | 192 |
| `group_algebra_keys` | 96 |
| `nonzero_group_algebra_keys` | 0 |
| `cancelling_pairs` | 96 |
| `bad_pairs` | 0 |
| `all_terms_abs_coeff` | 4 |
| `total_abs_before_cancellation` | 768 |
| `total_abs_after_cancellation` | 0 |

## Pattern

```text
signed_coefficient_patterns = {'(-4, 4)': 96}
source_destination_type_patterns = {"('split_pairs', 'split_pairs', (4, 4), (-1, 1))": 96}
```

## Boundary

Finite group-algebra cancellation certificate only; no full native B_n fingerprint theorem or classification theorem.

Prossima task: turn this finite cancellation into a short human-readable lemma, or return to final PDF/freeze routing.
