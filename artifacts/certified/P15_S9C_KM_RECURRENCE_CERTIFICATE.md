# P15-S9C K_m Recurrence Certificate

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S9C K_m recurrence certificate
Date: 2026-07-07
Status: PASS

## Decision

P15-S9C passes. The final odd pair-standard `k=1` residue

```text
K_m = (m-1)*lambda_prime(2m+1,1)
```

is positive for every `m>=2` as a displayed formula sequence. S9C itself does not identify this sequence with the true odd pair-standard `k=1` scalar; downstream S11A-K now supplies that closure.

This is still not a public theorem promotion by S9C alone; downstream S11A-K supplies the scalar-link closure, and S10 continues to govern public wording.

## Machine Artifact

- Script: `scripts/p15_s9c_km_recurrence_certificate.py`
- JSON: `artifacts/certified/P15_S9C_KM_RECURRENCE_CERTIFICATE.json`
- Schema: `p15.s9c.km_recurrence_certificate.v1`

Run command:

```text
python -B scripts/p15_s9c_km_recurrence_certificate.py --max-m 100 --write-json artifacts/certified/P15_S9C_KM_RECURRENCE_CERTIFICATE.json
```

## Certificate Shape

For `m>=4`, S9 gives

```text
lambda_prime(2m+1,1) = Phi_{2m}(B(t)^(m-4)*L_m(t)),
B(t)=2*t^2-4*t+1.
```

After the factorial-functional substitution `t=1/(2x)`, the script writes

```text
K_m = integral_0^infty exp(-x) * 2^(m+4) * A(x)^(m-4) * Q_m(x) dx,
A(x)=2*x^2-4*x+1.
```

It stores an explicit polynomial `S(m,x)` and verifies exactly over `Q[m,x]`:

```text
A*S_x + ((m-4)*A_x - A)*S
= A*sum_j p_j(m) 2^j A^j Q_{m+j}.
```

Boundary terms vanish because `S(m,0)=0` and `exp(-x)` kills infinity. This proves the order-3 recurrence for every `m>=4`; recurrence edge cases `m=2,3` are checked exactly.

## Checks Passed

```text
all_m_telescoping_certificate          PASS
edge_case_and_finite_sanity            PASS
all_m_K_m_positivity_from_recurrence   PASS
```

Key machine facts:

```text
R degrees: x=13, m=15, terms=217
S degrees: x=13, m=14, terms=192
polynomial residual terms: 0
S(m,0)=0: True
finite sanity m=2..100: no nonzero residuals
K_2..K_5 = 28, 14640, 5241648, 2587263104
```

The positivity induction uses `rho_m=K_{m+1}/K_m >= 4*m^2`. Base rows `m=2..8` pass. For `m>=7`, shifted-polynomial sign checks prove:

```text
p0(m)>0,
p1(m)<=0,
-p2(m)>0,
p3(m)>0,
D(m)>0.
```

Here

```text
D(m)=16*m^2*(m+1)^2*(-p2(m)) - p0(m)
     - 64*m^2*(m+1)^2*(m+2)^2*p3(m).
```

These inequalities imply the ratio step `rho_{m+2}>=4*(m+2)^2`, hence `K_m>0` for every `m>=2`.

## Closed At S9C

- The order-3 recurrence for `K_m` is certified by integration by parts.
- `K_m>0` is certified for all `m>=2`.
- Odd pair-standard `k=1` formula-sequence positivity is closed.
- Together with S9A/S9B and downstream S11A-K, the odd `V_pair-std` `k=1` true-scalar route is closed.
- Together with S4, S8, downstream S11A-H, and downstream S11A-K, the two ambient standard-channel scalar positivity chain is internally closed.

## Downstream Closed

- S11A-K formula-to-true-scalar closure is now closed downstream by `P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE`.

## Still Open

- Positive-only separator and full `B_n` fingerprint remain separate research tasks.
- `D~I` remains Type-D scout-only; no Type-D theorem is promoted.

Next integration task: manuscript polish and appendix routing for S8/S9C/S11A certificates.
