# P15-S4 Even-n Theorem Certificate

Project: Signed Reversal Rank Theorem (ID P15)
Gate: P15-S4 even-n theorem
Verdict: PASS
Date: 2026-07-06

## Gate requirement

Certify the even-n theorem for the two ambient standard channels:

```text
rank M_ref(B~I_n(k,k)) = n/2,
rank M_pair-std(B~I_n(k,k)) = n/2 - 1,
```

for even `n` and every `k>=1` with nonempty diagonal layer.

## Delivered artifacts

1. `docs/P15_S4_EVEN_N_THEOREM_2026_07_06.md` - theorem statement, proof, boundary, and replay linkage.
2. `scripts/p15_s4_even_theorem_certificate.py` - focused S4 certificate builder; imports the S3 exact engine and
   runs with `python -B` to avoid `__pycache__` writes.
3. `artifacts/certified/P15_S4_EVEN_N_THEOREM_CERTIFICATE.json` - machine-readable S4 certificate, UTF-8 no BOM.

## Observed command

```powershell
python -B scripts\p15_s4_even_theorem_certificate.py --write-json certified\P15_S4_EVEN_N_THEOREM_CERTIFICATE.json --limit-n 20
```

Observed result: `status: PASS`.

## Checks passed

| Check | Result | Range |
|---|---|---|
| Even exact channel replay | PASS: ranks, collapse, `a` closed form, trace/flip, and `lambda_pair` formula | full enumeration cases `n=4,6`, diagonal nonempty `k>=1` |
| Even closed-form positivity sweep | PASS: `a(n,k)>0` | even `n=4..20`, nonempty `k>=1` |
| Even permanent-domination replay | PASS: `E[#neg-fixed] <= 1` for every feasible pattern tested | all feasible patterns for even `n=4,6` |
| Proof checklist | PASS | S1/S2/S3 dependencies, collapse, permanent lemma, trace identities, boundary controls |

## Boundary carried forward

- S4 is an internal theorem certificate, not public promotion.
- The rank number and `P_+` projector remain classical matching-scheme bookkeeping.
- The P15 contribution remains collapse plus positivity for the signed-reversal layer.
- Odd `n`, `H_m`, odd `V_pair-std`, `C_2 x C_2` square, Type-D scout material, and full `B_n` fingerprint remain outside S4.

Verdict: PASS.  P15-S4 even-n theorem complete; next gate is P15-S5 odd-n structure.
