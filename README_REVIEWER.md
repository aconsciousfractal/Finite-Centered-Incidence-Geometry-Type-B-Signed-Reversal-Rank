# Reviewer Guide

Start with the claim boundary, then run the replay tiers you want to inspect.

Recommended reading order:

1. `docs/CLAIM_BOUNDARY.md`
2. `docs/CLAIM_LEDGER.md`
3. `paper/main.tex` and the split sources under `paper/sections/` and `paper/appendices/`
4. `REPRODUCE.md`
5. `docs/APPENDIX_BOUNDARY.md`

Replay tiers:

1. Core theorem replay: S4/S5/S7/S8/S9/S9C plus S11 scalar-link checks.
2. Type-D sublayer replay: `type_d/` G15/G16/G20/G22 (same two channels, same diagonal cells; `opus_g16` is the heavy one).
3. Separator add-on replay: S12 Gate A/B.
4. Optional full-fingerprint appendix replay: finite `n=8` zero mechanisms and selected fixed-plus audits.

The optional full-fingerprint appendix is intentionally not part of the main theorem claim.

Note on dated certificates. Some dated per-gate documents and frozen certificates under `docs/` and
`artifacts/certified/` predate the type-D parity-sublayer import and still carry per-gate language such as
"`D~I` remains scout-only; no Type-D theorem is promoted." These are accurate *per gate* — those gates
(S4/S5/S6/S7/S8/S9/S9C/S11) do not themselves promote the type-D theorem — and are kept as historical records.
The current live claim boundary is `docs/CLAIM_BOUNDARY.md`, `PUBLIC_CLAIM_BOUNDARY.md`, and the type-D sublayer
material in `paper/sections/11_type_d_parity_sublayer.tex` and `type_d/`.
