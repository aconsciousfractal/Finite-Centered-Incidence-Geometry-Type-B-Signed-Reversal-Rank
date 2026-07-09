# P15-S10 Public Theorem Statement

Date: 2026-07-07
Status: PUBLIC-SAFE AFTER S10 CHECK; POST-S12 AND POST-S14 ADDENDA RECORDED BELOW

## Scoped Theorem

Let `B_n` be the signed permutation group and let `rho(i)=n-1-i`. Define

```text
B~I_n(k,l) = { (pi,eps) in B_n : #{i: pi(i)=i, eps_i=+}=k,
                                 #{i: pi(i)=rho(i), eps_i=+}=l }.
```

For every `n>=2` and every `k>=1` such that the diagonal layer `B~I_n(k,k)` is nonempty, the two ambient standard-channel aggregate operators satisfy

```text
rank M_ref(B~I_n(k,k)) = ceil(n/2)
rank M_pair-std(B~I_n(k,k)) = ceil(n/2)-1.
```

This is the public-safe rank theorem statement after S10.

## Post-S12 Separator Add-On

The following additional statement is public-safe after S12 Gate A/B:

For `X_n=B~I_n(1,1)`, in the positive-only signed-reversal layer, and for every `n>=4`, `X_n` is neither a Reiner-style signed-poset linear-extension set nor an exact left-translate tiler of `B_n`. The case `n=3` is recorded separately as an exact left-translate tiling exception.

This add-on is scoped only to `X_n=B~I_n(1,1)` and does not claim a full Fourier fingerprint, a classification theorem, or a positive tiling theorem.

## Post-S14 Type-D Add-On

Date: 2026-07-09. The following additional statement is public-safe after the P15-D extension import (S14):

For every `n>=2` and every `k>=1` such that the type-D parity sublayer `D~I_n(k,k)=B~I_n(k,k) ∩ D_n` is nonempty (`D_n=ker chi`, `chi(pi,eps)=prod_i eps_i` the even-sign subgroup), the same two ambient standard-channel aggregate operators satisfy

```text
rank M_ref(D~I_n(k,k)) = ceil(n/2)
rank M_pair-std(D~I_n(k,k)) = ceil(n/2)-1.
```

A sign-character split reduces the parity restriction to a sign-weighted correction scalar, closed by a factorized generating function (even and odd reference) and by a direct type-D trace identity (odd pair-standard); replay in `type_d/`. This add-on is confined to the same two ambient standard channels and the same diagonal cells; it claims no other reflection subgroup, channel, or general Weyl statement.

## Public-Safe Novelty Framing

The contribution is the object-specific collapse of the signed-reversal diagonal layer onto the matching-scheme idempotent, together with the scalar positivity proofs that make the ranks exact.

Enumeration is classical background. The rank number, reversal projector, and matching-scheme multiplicities are classical bookkeeping and must be cited rather than claimed as new.

## Not Claimed

- Beyond the Post-S14 type-D parity-sublayer add-on above (same two channels, same diagonal cells), no theorem about other `D~I` cells or off-channel type-D layers.
- No general Weyl-group theorem.
- No affine, positive-tiling, classification, or full fingerprint theorem; no separator theorem beyond the scoped S12 statement for `X_n=B~I_n(1,1)`.
- No novelty claim for enumeration, rank numbers, projectors, or matching-scheme multiplicities.
