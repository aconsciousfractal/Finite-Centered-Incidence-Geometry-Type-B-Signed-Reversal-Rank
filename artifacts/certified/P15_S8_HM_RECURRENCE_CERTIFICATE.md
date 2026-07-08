# P15-S8 H_m Recurrence Certificate

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S8 H_m recurrence certificate
Date: 2026-07-07
Status: PASS

## Decision

P15-S8 passes for the displayed `H_m` formula sequence. The order-3 recurrence for that sequence is certified for all `m>=2`, and the elementary ratio induction proves `H_m>0` for all `m>=2`.

This proves positivity of the displayed formula sequence used in the odd `V_ref`, `k=1` route:

```text
a(2m+1,1) = a(2m,1) + H_m,
H_m > 0,
a(2m,1) > 0.
```

S8 does not by itself prove the S11A-H formula-to-true-scalar link and does not promote a public theorem.

## Machine Artifact

- Script: `scripts/p15_s8_hm_recurrence_certificate.py`
- JSON: `artifacts/certified/P15_S8_HM_RECURRENCE_CERTIFICATE.json`
- Schema: `p15.s8.hm_recurrence_certificate.v1`

Run command:

```text
python -B scripts/p15_s8_hm_recurrence_certificate.py --max-m 100 --write-json artifacts/certified/P15_S8_HM_RECURRENCE_CERTIFICATE.json
```

## Certificate Shape

For `m>=3`, write the factorial functional as an integral:

```text
H_m = integral_0^infty exp(-x) * 2^m * B(x)^(m-3) * Cbar_m(x) dx,
B(x)=2x^2-4x+1.
```

The script stores an explicit polynomial `S(m,x)` and verifies exactly over `Q[m,x]`:

```text
B*S_x + ((m-3)*B_x - B)*S = B*R,
R = sum_{j=0}^3 p_j(m) 2^j B^j Cbar_{m+j}.
```

Therefore

```text
exp(-x)2^m B^(m-3)R = d/dx(exp(-x)2^m B^(m-3)S).
```

Boundary terms vanish: `S(m,0)=0`, and `exp(-x)` kills the infinity endpoint. This proves the recurrence for every `m>=3`. The `m=2` edge case is checked exactly.

## Checks Passed

```text
all_m_telescoping_certificate                  PASS
edge_case_and_finite_sanity                    PASS
all_m_H_m_positivity_from_recurrence           PASS
```

Key machine facts:

```text
R degrees: x=11, m=10, terms=125
S degrees: x=11, m=9, terms=107
polynomial residual terms: 0
S(m,0)=0: True
m=2 recurrence residual: 0
finite sanity m=2..100: no nonzero residuals
```

The positivity induction checks signs of `p0`, `p1`, `-p2`, `p3`, positivity of `D(m)`, and base ratios for `m=2,3,4`.

## Closed At S8

- The order-3 recurrence for `H_m` is no longer just guessed/held-out; it has a checked integration-by-parts certificate.
- `H_m>0` is certified for all `m>=2`.
- S8 itself does not identify this displayed formula sequence with the true odd `V_ref`, `k=1` scalar residual; downstream S11A-H now supplies that closure.

## Downstream Closed

- S11A-H formula-to-true-scalar closure is now closed downstream by the S11A-H edge-insertion verifier.

## Still Open

- Odd pair-standard downstream closure and S11A-K formula-to-true-scalar closure are now closed downstream by S9/S9C/S11A-K.
- Positive-only separator hunt [subsequently closed by the S12 separator theorem] and the full native `B_n` fingerprint [remains open by design].
- `D~I` remains Type-D scout-only; no Type-D theorem is promoted.

Verdict: PASS. Next gate: P15-S9 manuscript shape, odd pair-standard positivity, and separator/fingerprint boundary work.
