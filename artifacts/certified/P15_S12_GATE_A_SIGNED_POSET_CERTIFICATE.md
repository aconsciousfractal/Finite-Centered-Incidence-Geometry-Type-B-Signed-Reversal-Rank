# P15 S12 Gate A Signed-Poset Certificate

Schema: `p15.s12_gate_a_signed_poset_certificate.v1`
Status: **GATE_A_CERTIFICATE**

Convention: `T(w)=-w_n<...<-w_1<0<w_1<...<w_n`.

## Summary

- `n4_ordered_nonzero_pairs`: `56`
- `n4_all_nonzero_precedences_killed`: `True`
- `n4_ordered_pairs_with_zero`: `72`
- `n4_all_precedences_with_zero_killed`: `True`
- `n4_X_size`: `40`
- `n4_distinct_witnesses_with_zero`: `3`
- `sign_extension_all_specs_n5_n7`: `True`
- `sign_extension_all_one_label_specs_n4_n7`: `True`
- `sign_extension_n4_failure_count`: `2`
- `gate_a_status`: `PASS_FOR_SIGNED_POSET_ROUTE`

## n=4 Finite Witnesses

For every ordered signed pair `(a,b)`, including pairs with `0`, the certificate stores a witness `w in X_4` with `b` before `a` in `T(w)`. This kills the possible common precedence `a<b`.

| a | b | witness w |
|---:|---:|---|
| -4 | -3 | `[1, -3, 2, -4]` |
| -4 | -2 | `[1, -3, 2, -4]` |
| -4 | -1 | `[1, -3, 2, -4]` |
| -4 | 0 | `[1, -3, 2, -4]` |
| -4 | 1 | `[1, -3, 2, -4]` |
| -4 | 2 | `[1, -3, 2, -4]` |
| -4 | 3 | `[1, -3, 2, -4]` |
| -4 | 4 | `[1, -3, 2, -4]` |
| -3 | -4 | `[-1, -3, 2, 4]` |
| -3 | -2 | `[-1, -3, 2, 4]` |
| -3 | -1 | `[-1, -3, 2, 4]` |
| -3 | 0 | `[-1, -3, 2, 4]` |
| -3 | 1 | `[-1, -3, 2, 4]` |
| -3 | 2 | `[-1, 3, -2, 4]` |
| -3 | 3 | `[-1, -3, 2, 4]` |
| -3 | 4 | `[1, -3, 2, -4]` |
| -2 | -4 | `[-1, -3, 2, 4]` |
| -2 | -3 | `[-1, 3, -2, 4]` |
| -2 | -1 | `[-1, 3, -2, 4]` |
| -2 | 0 | `[-1, 3, -2, 4]` |
| -2 | 1 | `[-1, 3, -2, 4]` |
| -2 | 2 | `[-1, 3, -2, 4]` |
| -2 | 3 | `[-1, 3, -2, 4]` |
| -2 | 4 | `[1, -3, 2, -4]` |

Only the first 24 witness rows are printed here; the JSON contains all rows.

## Sign-Extension Smoke Check

| n | total sign specs | failures | one-label failures | all extend? | failure sample |
|---:|---:|---:|---:|---|---|
| 4 | 32 | 2 | 0 | False | `[[-1, -4], [-2, -3]]` |
| 5 | 50 | 0 | 0 | True | `[]` |
| 6 | 72 | 0 | 0 | True | `[]` |
| 7 | 98 | 0 | 0 | True | `[]` |

## Interpretation

Gate A is closed for the signed-poset route: `n=4` is finite-certified, one-label sign freedom kills zero-relations in `n=4..7`, and the `n>=5` proof uses the sign-extension lemma written in the companion proof note.

## Reproducibility

Run from the project root:

`python -B scripts/p15_s12_gate_a_signed_poset_certificate.py --write-json certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.json --write-md certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.md`
