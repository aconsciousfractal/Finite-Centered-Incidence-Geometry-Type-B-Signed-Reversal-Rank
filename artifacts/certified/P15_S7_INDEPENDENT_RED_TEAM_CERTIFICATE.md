# P15-S7 Independent Red-Team Certificate

Project: Signed Reversal Rank Theorem
ID: P15
Gate: P15-S7 independent red-team recompute
Date: 2026-07-06
Status: PASS

## Decision

P15-S7 passes. A second verifier, `scripts/p15_s7_independent_red_team.py`, recomputes the active S3-S6 surface without importing `p15_s3_exact_engine.py` or any S4-S6 certificate script.

S7 is not a new infinite theorem. It is a red-team gate: it asks whether an independent implementation can reproduce the finite exact channel ranks, count identities, and boundary discipline already used by S3-S6. It can.

## Machine Artifact

- Script: `scripts/p15_s7_independent_red_team.py`
- JSON: `artifacts/certified/P15_S7_INDEPENDENT_RED_TEAM_CERTIFICATE.json`
- Schema: `p15.s7.independent_red_team.v1`

Run command:

```text
python -B scripts/p15_s7_independent_red_team.py --write-json artifacts/certified/P15_S7_INDEPENDENT_RED_TEAM_CERTIFICATE.json
```

## Checks Passed

1. Independent channel recompute: raw signed permutations rebuild `B~I_n(k,k)`, `M_ref`, and `M_pair` for `n=3..7`, all nonempty `k>=1` diagonal rows.
2. Exact rank check: `rank M_ref=ceil(n/2)` and `rank M_pair-std=ceil(n/2)-1` on all checked rows.
3. Collapse check: even rows satisfy `M_ref=2a.P_+`; odd rows satisfy `M_ref=2a.P_+^{pairs}+a_c.E_center`.
4. Base marker: `n=3,k=1` is rank-correct with `a=-1`, `a_c=2`; S7 does not reinterpret that as positivity.
5. Enumeration boundary: A000354/type-B derangement prefix `[1,1,5,29,233,2329,27949,391285,6260561]` is reproduced, and `m_n(k)=C(n,k)A000354(n-k)` is brute-forced for `n=0..7`.
6. Square/count red-team: deterministic integer evaluations verify `F_square=per(W)` for `n=2..7`; the `q=-1` collapse is checked at integer points; exact row checks verify old `BI` and new `B~I` marginals, old `DI` recovery, and even `D~I` correction for `n<=6`.
7. Boundary scan: active governance/docs are clean of the stale if-and-only-if `H_m` wording and still state the A000354/A007016, matching-scheme, no-public-claim, and Type-D scout boundaries.

## Boundary Correction Made Before Final S7

S7 found one stale S1 wording in `docs/P15_S1_OBJECT_ADMISSION_AND_CONVENTIONS_2026_07_06.md`. The old O6 line used an if-and-only-if phrasing that was too strong for the current gate.

The active wording now records the S5 reduction correctly:

```text
a(2m+1,1) = a(2m,1) + H_m
```

and leaves `H_m>0` as the S8 certificate obligation.

## Still Open

- P15-S8: certify the `H_m` recurrence / creative-telescoping step for odd `k=1`.
- General odd `V_pair-std` positivity for `lam_prime` and `mu`.
- P15-S9: manuscript/public-boundary shaping and separator/full-fingerprint scouting.
- `D~I` remains Type-D scout-only; no Type-D theorem is promoted.

Verdict: PASS. Next gate: P15-S8 `H_m` recurrence certificate.
