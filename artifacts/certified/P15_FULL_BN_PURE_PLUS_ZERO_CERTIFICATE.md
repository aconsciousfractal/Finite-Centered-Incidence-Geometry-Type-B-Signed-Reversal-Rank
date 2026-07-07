# P15 Full Bn Pure-Plus Zero Block Certificate

Status: **PASS**

## Claim

```text
M_{(2.2.2.1.1 | empty)}(B~I_8(3,3)) = 0 in the rational S_8 Specht model
```

## Interpretation

This is the remaining `ZERO` row in the n=8, k=3 full native `B_n` atlas after the `(1.1|nu)`, `nu |- 6`, zero family is removed.
Because `lambda_minus` is empty, signs carry no nontrivial character. The coefficient of an underlying permutation `pi in S_8` is exactly the number of sign choices producing three positive identity hits and three positive reversal hits.
Only 102 underlying permutations have nonzero coefficient. Summing their exact Young-seminormal `S_8` matrices on shape `(2,2,2,1,1)`, with these integer weights, gives the zero rational matrix.

## Exact Audit

| field | value |
|---|---:|
| `specht_dim` | 28 |
| `nonzero_underlying_permutations` | 102 |
| `layer_size` | 480 |
| `common_denominator_seen` | 2160 |
| `max_abs_numerator_seen` | 245 |
| `nonzero_entries` | 0 |

## Weight Patterns

```text
{'(3, 3, 4)': 96, '(4, 4, 16)': 6}
```

## Boundary

Finite exact pure-plus zero-block witness only; no full native B_n fingerprint theorem or classification theorem.

Prossima task: record the two zero mechanisms together in the n=8 defect note, then decide whether to search for a conceptual proof of the pure-plus zero or keep it as a finite certified witness.
