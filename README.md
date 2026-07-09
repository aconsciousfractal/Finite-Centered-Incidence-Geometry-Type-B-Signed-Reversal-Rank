# Exact Standard-Channel Rank of a Signed-Reversal Incidence Layer

This repository contains the public source, certificates, and replay scripts for the P15 signed-reversal rank theorem.

## Main Claim

For the signed-reversal layer `B~I_n(k,k)` in the signed permutation group `B_n`, the two ambient standard-channel aggregate operators have ranks

```text
rank M_ref(B~I_n(k,k))      = ceil(n/2)
rank M_pair-std(B~I_n(k,k)) = ceil(n/2)-1
```

for every nonempty diagonal row with `n>=2` and `k>=1`.

The rank number and reversal projector are classical matching-scheme bookkeeping. The contribution is the object-specific collapse and scalar positivity/reduction for this noncentral signed-reversal layer.

## Type-D Parity Sublayer

For the same diagonal rows, the type-D parity sublayer `D~I_n(k,k) = B~I_n(k,k) ∩ D_n`, where `D_n = ker(chi)` is the even-sign subgroup of `B_n` and `chi(pi,eps) = prod_i eps_i`, has the same two standard-channel ranks:

```text
rank M_ref(D~I_n(k,k))      = ceil(n/2)
rank M_pair-std(D~I_n(k,k)) = ceil(n/2)-1
```

A sign-character split writes the type-D aggregate as the average of the sign-blind and the sign-weighted aggregates; the sign-weighted correction is controlled by an explicit factorized generating function and, in the odd pair-standard case, by a direct type-D trace identity. Replay scripts and certificates are in `type_d/`. The result is confined to the same two ambient standard channels and the same diagonal cells as the full-layer theorem.

## Additional Scoped Result

For `X_n=B~I_n(1,1)` and `n>=4`, the positive-only layer is neither a Reiner-style signed-poset linear-extension set nor an exact left-translate tiler of `B_n`. The `n=3` row is recorded as an exact-tiling exception.

## Appendix Boundary

The finite native `B_n` fingerprint appendix records `n=8` evidence only. It explains the twelve `n=8,k=3` zero rows and selected finite extra-defect families, but it is not a full native `B_n` fingerprint theorem or classification.

Current finite atlas accounting:

```text
114 defect rows total
44 accounted rows
70 open rows
```

## Repository Layout

```text
paper/main.tex                 manuscript driver
paper/Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer.pdf  titled reviewer PDF
paper/abstract.tex             abstract text
paper/sections/                main-text sections
paper/appendices/              reproducibility and finite-atlas appendices
paper/refs.bib                 bibliography database
scripts/                       deterministic replay scripts
type_d/                        type-D parity sublayer replay scripts and certificates
artifacts/certified/           certified JSON/Markdown artifacts
artifacts/external_n8_atlas/   finite n=8 atlas inputs
docs/                          claim boundary and source-lock docs
REPRODUCE.md                   replay and paper-build commands
README_REVIEWER.md             suggested reviewer path
```

## Quick Smoke Test

```text
pip install -r requirements.txt
python -B scripts/p15_s4_even_theorem_certificate.py --write-json artifacts/certified/P15_S4_EVEN_N_THEOREM_CERTIFICATE.json
```

This command checks only the even-row core certificate. For the full theorem replay, use `REPRODUCE.md`, especially the Core Theorem Replay and Scalar-Link Replay sections. See `REPRODUCE.md` also for the `pdflatex`/`bibtex` paper build.
