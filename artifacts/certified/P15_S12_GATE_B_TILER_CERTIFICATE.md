# P15 S12 Gate B Tiler Certificate

Schema: `p15.s12_gate_b_tiler_certificate.v1`
Status: **GATE_B_CERTIFICATE_CLOSED**

## Summary

- `finite_n3_is_divisible_exception`: `True`
- `n3_exact_left_translate_tiling_exists`: `True`
- `n3_exact_left_translate_tile_count`: `24`
- `finite_windows_ok_n4_n21`: `True`
- `finite_divisibility_fails_n4_n21`: `True`
- `symbolic_truncation_all_positive`: `True`
- `tail_probe_signs_ok_to_max`: `True`
- `tail_probe_decreasing_to_max`: `True`
- `max_tail_probe_n`: `120`
- `gate_b_status`: `PASS: finite windows, symbolic truncations, and tail monotonicity proof closed`

## Finite Window Table

| n | floor(G/N) | G mod N | window | ok? | divides? |
|---:|---:|---:|---|---|---|
| 3 | 24 | 0 | n=3 exception | False | True |
| 4 | 9 | 24 | 9 < G/N < 10 (finite_n4) | True | False |
| 5 | 10 | 200 | 10 < G/N < 11 (odd) | True | False |
| 6 | 11 | 1200 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 7 | 10 | 43600 | 10 < G/N < 11 (odd) | True | False |
| 8 | 11 | 173760 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 9 | 10 | 12293600 | 10 < G/N < 11 (odd) | True | False |
| 10 | 11 | 45031040 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 11 | 10 | 5462495680 | 10 < G/N < 11 (odd) | True | False |
| 12 | 11 | 16760726400 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 13 | 10 | 3459331026560 | 10 < G/N < 11 (odd) | True | False |
| 14 | 11 | 8351814831360 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 15 | 10 | 2948739784433920 | 10 < G/N < 11 (odd) | True | False |
| 16 | 11 | 5160953021593600 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 17 | 10 | 3250887638626895360 | 10 < G/N < 11 (odd) | True | False |
| 18 | 11 | 3546093092226600960 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 19 | 10 | 4498934950931210030080 | 10 < G/N < 11 (odd) | True | False |
| 20 | 11 | 1978435270261426268160 | 11 < G/N < 12 (even_6_to_20) | True | False |
| 21 | 10 | 7634524049443036808284160 | 10 < G/N < 11 (odd) | True | False |

## n=3 Exact Tiling Witness

- `left_translate_tiling_exists`: `True`
- `target_tile_count`: `24`
- `chosen_left_translate_reps`: `[[-2, 1, 3], [1, 3, 2], [-2, 1, -3], [-1, 3, 2], [-1, 2, -3], [1, 2, -3], [-1, 3, -2], [2, 1, 3], [2, 1, -3], [-1, -3, -2], [-2, -1, 3], [-1, -2, 3], [-1, 2, 3], [-1, -2, -3], [-2, -1, -3], [2, -1, -3], [1, -2, -3], [1, 3, -2], [1, -3, 2], [1, -2, 3], [1, -3, -2], [1, 2, 3], [2, -1, 3], [-1, -3, 2]]`

## Symbolic Truncation Certificates

| certificate | shift | degree | all coeffs positive? | coefficients descending |
|---|---:|---:|---|---|
| `even_m_ge_11_lower_L9_minus_1_over_11` | 11 | 4 | True | `[29952, 753728, 6383328, 18796492, 5119335]` |
| `even_m_ge_11_upper_1_over_10_minus_U6` | 11 | 3 | True | `[48, 1528, 15693, 52493]` |
| `odd_m_ge_11_lower_L9_minus_1_over_11` | 11 | 6 | True | `[479232, 32334080, 862589824, 11849156608, 89258599064, 351693457047, 568489477815]` |
| `odd_m_ge_11_upper_1_over_10_minus_U6` | 11 | 4 | True | `[48, 1960, 29845, 201040, 505827]` |

## Odd Tail Coefficient Reduction

`S_{2m+1}(-x)=x f(x)^(m-2) Q_m(x)`
`Q_m` coefficients x0..x4: `['-1', '4*(m - 2)*(m + 1)', '4*(3*m**2 - 3*m - 5)', '4*(3*m**2 - 3*m - 4)', '4*(m**2 - m - 1)']`
`2((2m+1)C_m-xC_m')-C_m/x=f(x)^(m-3)P_m(x)`
`P_m` coefficients x0..x6: `['1', '-4*(m**2 - 3)', '2*(8*m**3 - 26*m**2 - 2*m + 27)', '4*(5*m - 4)*(4*m**2 - 5*m - 7)', '4*(36*m**3 - 61*m**2 - 31*m + 27)', '8*(14*m**3 - 23*m**2 - 9*m + 6)', '8*(4*m**3 - 7*m**2 - m + 1)']`
Tail lemma status: `closed in P15_S12_GATE_B_TILER_PROOF_2026_07_07.md via Newton/log-concavity`

## Interpretation

This certificate records the finite range, the n=3 exact tiling witness, symbolic truncation inequalities, and algebraic odd-tail reduction used by the closed Gate B proof note.

## Reproducibility

Run from the project root:

`python -B scripts/p15_s12_gate_b_tiler_certificate.py --max-tail-probe-n 120 --write-json artifacts/certified/P15_S12_GATE_B_TILER_CERTIFICATE.json --write-md artifacts/certified/P15_S12_GATE_B_TILER_CERTIFICATE.md`
