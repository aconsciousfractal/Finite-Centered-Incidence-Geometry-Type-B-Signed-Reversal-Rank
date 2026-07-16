# Type-D Parity Sublayer Certificates

Deterministic exact-arithmetic replay scripts and JSON certificates for the type-D
parity sublayer `D~I_n(k,k) = B~I_n(k,k) ∩ D_n` (Section "The Type-D Parity Sublayer"
of the manuscript, Theorem `thm:typed`). `D_n = ker(chi)` is the even-sign subgroup of
`B_n`, `chi(pi,eps) = prod_i eps_i`.

These scripts use only the Python standard library (`fractions`, `math`, `itertools`);
no NumPy/SymPy is required. Each accepts `--no-write` (skip writing the JSON) and
`--out <path>` (override the output path), and each reports a top-level `status` field
(`PASS`/`FAIL`) plus `pass`, `checks_passed`, `checks_total`. Exit code is `0` on PASS.

| Script | Certifies | JSON artifact |
|---|---|---|
| `p15d_g15_C_closure_verify.py` | odd sign-localization: `s^chi >= 0` on the band `3k > 2m` (sign-definiteness, per-index ratio bound, single-orphan structure) | `P15D_G15_C_CLOSURE.json` |
| `p15d_g20_oatail_symbolic_verify.py` | orphan-tail closed forms `R_A, R_B` and the reduction of `s^chi > 0` to two nonnegative quartics | `P15D_G20_OATAIL_SYMBOLIC.json` |
| `p15d_g22_psdom_typeD_trace_verify.py` | odd pair-standard `lambda'^D > 0` via the direct type-D trace identity (finite replay `n=5,7`; `k=1` via the `K_m > D_m` comparison) | `P15D_G22_PSDOM_TYPED_TRACE.json` |
| `p15d_g16_pairstd_closure_verify.py` | closed form of the sign-weighted scalar `lambda'^chi` and the auxiliary domination cross-checks (`n = 5,7,9,11,13`) | `P15D_G16_PAIRSTD_CLOSURE.json` |

## Dependency

`p15d_g22_psdom_typeD_trace_verify.py` imports the shared exact engine
`p15_s3_exact_engine.py` from the repository's top-level `scripts/` directory
(`../scripts/`); the other three scripts are self-contained. Run from anywhere; the
engine path is resolved relative to the script location.

## Replay

```text
python -B type_d/p15d_g15_C_closure_verify.py  --out type_d/P15D_G15_C_CLOSURE.json
python -B type_d/p15d_g20_oatail_symbolic_verify.py --out type_d/P15D_G20_OATAIL_SYMBOLIC.json
python -B type_d/p15d_g22_psdom_typeD_trace_verify.py --out type_d/P15D_G22_PSDOM_TYPED_TRACE.json
python -B type_d/p15d_g16_pairstd_closure_verify.py --out type_d/P15D_G16_PAIRSTD_CLOSURE.json
```

`p15d_g16_pairstd_closure_verify.py` recomputes exact permanents up to `n = 13` and is
the heavy one (several minutes); the other three complete in well under a minute.

## Boundary

These certificates back the type-D rank equalities on the same two ambient standard
channels and the same diagonal cells as the full-layer theorem. They do not add a full
native `B_n` fingerprint theorem, a classification theorem, or any general Weyl claim.
