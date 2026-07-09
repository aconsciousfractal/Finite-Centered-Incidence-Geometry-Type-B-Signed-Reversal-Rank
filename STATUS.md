# Public Boundary Summary

This public repository package promotes the scoped two-standard-channel rank theorem, its type-D parity sublayer, and the scoped S12 separator theorem.
The type-D parity sublayer `D~I_n(k,k) = B~I_n(k,k) ∩ D_n` carries the same two standard-channel ranks, on the same channels and diagonal cells.
There is no public theorem claim for general Weyl, affine, or full native B_n fingerprint classification.
The signed enumeration is classical background and is not A007016 for the signed enumeration.
The rank number and reversal projector are matching-scheme bookkeeping.

## Closed Theorem Chain

The standard-channel theorem is internally closed by the following replay/proof chain:

- S4 even theorem: matching-scheme collapse and scalar positivity for even `n`.
- S5 odd reductions: odd rows reduce to even inputs plus two `k=1` residuals.
- S8 plus S11A-H: the `H_m` recurrence and the true-scalar link close the odd reference residual.
- S9/S9C plus S11A-K: the `K_m` recurrence and true-scalar/product-Phi links close the odd pair-standard residual.
- S11D smoke replay: independent exact rank checks on small rows confirm the end-to-end rank numbers.
- Type-D parity sublayer: a sign-character split (`type_d/`) reduces the parity restriction to a sign-weighted correction scalar, closed by a factorized generating function (even and odd reference) and a direct type-D trace identity (odd pair-standard); certificates `P15D_G15/G16/G20/G22`.

## Add-Ons And Boundaries

The S12 positive-only separator theorem is a separate scoped add-on for `X_n=B~I_n(1,1)`.
The finite native `B_n` appendix is boundary evidence only; it is not a full fingerprint theorem or classification.
The type-D parity sublayer `D~I` is a proved rank result on the same two standard channels and diagonal cells (`type_d/`); it is not a general Weyl or affine claim.

Current status: public repository package; replay tiers pass locally; split paper PDF generated for human review.
