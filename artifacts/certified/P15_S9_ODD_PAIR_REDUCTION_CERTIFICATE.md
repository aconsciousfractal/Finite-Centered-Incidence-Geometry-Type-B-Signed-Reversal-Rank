# P15-S9 Odd Pair-Standard Reduction Certificate

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S9 odd pair-standard reduction
Date: 2026-07-07
Status: PASS

## Decision

P15-S9A/S9B pass. In the odd pair-standard channel, `mu>0` is closed for every nonempty odd diagonal row, and `lambda_prime>0` is closed for every nonempty odd diagonal row with `k>=2`.

S9 itself isolates the last odd pair-standard residue

```text
K_m = (m-1) lambda_prime(2m+1,1) > 0,  m>=2.
```

That residue is now closed downstream by the P15-S9C K_m recurrence certificate. The current certified chain has no remaining odd pair-standard positivity residue.

## Machine Artifact

- Script: `scripts/p15_s9_odd_pair_reduction_certificate.py`
- JSON: `artifacts/certified/P15_S9_ODD_PAIR_REDUCTION_CERTIFICATE.json`
- Schema: `p15.s9.odd_pair_reduction_certificate.v1`

Run command:

```text
python -B scripts/p15_s9_odd_pair_reduction_certificate.py --max-m 10 --max-formula-m 10 --write-json artifacts/certified/P15_S9_ODD_PAIR_REDUCTION_CERTIFICATE.json
```

## Certificate Shape

For `n=2m+1`, `N=2m`, write

```text
d_p = |B~I_N(p,p)|,
e_k = |B~I_N(k+1,k)|.
```

The center split gives

```text
|X| = d_{k-1}+d_k + 2*((N-2k)*d_k + 2*(k+1)*e_k)
g_center = d_{k-1}+d_k.
```

Therefore

```text
mu = ((N+1)*g_center - |X|)/N
   = (d_{k-1}-d_k) + (4/N)*(k*d_k-(k+1)*e_k)
   = a_c + 4*a_even.
```

S5 supplies `a_c>0` on nonempty odd rows, and S4 supplies `a_even>0` on the even rows where it is needed. The top odd rows are covered by `a_c`. Hence `mu>0`.

The trace identity is

```text
(m-1)*lambda_prime = (k-1)*|X| + NegFixSum - mu.
```

Since `g_center<=|X|`, one has `mu<=|X|`. Thus for `k>=3`, the right side is at least `(k-2)*|X|+NegFixSum>0`. For `k=2`, the certificate records an explicit signed permutation witness with one negative fixed point for every odd `n>=5`; hence `NegFixSum>0` and the right side is positive.

## Checks Passed

```text
mu_bridge_symbolic_identity      PASS
mu_bridge_closed_form_sweep      PASS   odd n=3..21
lambda_trace_finite_replay       PASS   exact n=3,5,7
lambda_k_ge_2_closure            PASS   witness odd n=5..21
lambda_k1_residue_target         PASS   formula values m=2..10
```

Key machine facts:

```text
mu symbolic residual terms: 0
lambda_prime(5,1) = 28,        K_2 = 28
lambda_prime(7,1) = 7320,      K_3 = 14640
lambda_prime(9,1) = 1747216,   K_4 = 5241648
lambda_prime(11,1)= 646815776, K_5 = 2587263104
```

For `m>=4`, the script checks that the bivariate coefficient formula and the univariate residue formula `Phi_{2m}(B^(m-4)L_m)` agree through `m=10`. P15-S9C then certifies the resulting `K_m` sequence for all `m>=2`.

## Closed At S9

- Odd pair-standard `mu` positivity is closed by the bridge `mu=a_c+4*a_even` and the S4/S5 positivity inputs.
- Odd pair-standard `lambda_prime` positivity is closed for every nonempty row with `k>=2`.
- The pair-standard odd problem is reduced to one scalar sequence, `K_m`.

## Downstream Closure

- P15-S9C certifies `K_m>0` for all `m>=2` by an S8-style recurrence and integration-by-parts certificate plus ratio induction.
- Therefore the current certified chain closes the odd pair-standard positivity problem.
- No Type-D, full `B_n` fingerprint, positive-only separator, or public theorem claim is promoted by S9.

Verdict: PASS for S9A/S9B reduction. Downstream closure: P15-S9C `K_m` positivity certificate PASS.
