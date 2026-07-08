# P15 Public Repository Audit

Status: **PUBLIC REPOSITORY MATERIALIZED / SMOKE REPLAY PASS / RED-TEAM POLISH PASS / TITLED PDF GENERATED / HEAVY OPTIONAL TIER NOT RUN**

Date: 2026-07-07

## Repository Root

```text
<repository root>
```

## Structure

The repository follows the Type-A public repo shape:

```text
README.md
README_REVIEWER.md
REPRODUCE.md
requirements.txt
CITATION.cff
LICENSE
SHA256SUMS.txt
paper/main.tex
docs/
scripts/
artifacts/certified/
artifacts/external_n8_atlas/
experiments/external_agents_2026_07_07/
figures/
```

Scripts are intentionally kept flat in `scripts/` because several optional full-fingerprint scripts import helper modules by filename.

## Replay Executed From Repository Package

Core theorem and scalar-link smoke replay:

```text
S4 even theorem: PASS
S5 odd structure: PASS
S7 independent red-team: PASS after adding public boundary stubs expected by the script
S8 H_m recurrence: PASS
S9 odd pair reduction: PASS
S9C K_m recurrence: PASS
S11A-H true scalar: PASS
S11A transfer: PASS
S11A-K local: PASS
S11A-K global: PASS
S11 raw rank smoke: PASS
```

Separator replay:

```text
S12 Gate A signed-poset: PASS
S12 Gate B tiler: PASS
```

Optional full-fingerprint smoke replay:

```text
(1.1|nu) group-algebra zero family: PASS
pure-plus Johnson degree-drop: PASS
pure-plus direct rational Specht zero: PASS
lambda_plus=(4,1), n=8,k=3: PASS
```

Paper build:

```text
pdflatex paper/main.tex: PASS on two passes, output in C:\tmp\p15_public_repo_latex
```

Static checks:

```text
No private absolute host paths found in repository content.
No __pycache__ directories found.
No LaTeX build byproducts found inside the repository package.
Python AST parse of all scripts: PASS.
SHA256SUMS regenerated after replay and script cleanup.
```

## Not Run Automatically

The exact rational `T_1/T_2/T_4` direct-rank audits were not run in the smoke pass. They are intentionally listed as a heavier optional tier in `REPRODUCE.md`:

```text
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T1 ...
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T2 ...
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T4 ...
```

## Fixes Made During Packaging

1. Added root public boundary stubs expected by `p15_s7_independent_red_team.py`:
   `STATUS.md`, `ROUTE_PLAN.md`, `CLAIM_LEDGER.md`, `PROOF_OBLIGATIONS.md`, `PUBLIC_CLAIM_BOUNDARY.md`, and selected gate stub docs.
2. Corrected `REPRODUCE.md`: `p15_s11a_scalar_link_symbolic_verifier.py` and `p15_s11a_k_true_scalar_bivariate_verifier.py` accept `--write-json` only, not `--write-md`.
3. Sanitized stale next-task strings in optional full-fingerprint scripts copied into the repository package.
4. Removed internal packaging notes from the public repository content.


## Postscript: Paper First Pass

After the public-package smoke audit, `paper/main.tex` was synced from the first reader-facing manuscript pass. The revised repository paper compiles cleanly to `C:\tmp\p15_public_repo_latex\Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer.pdf` in 11 pages, and `SHA256SUMS.txt` was regenerated after the sync.

## Postscript: Split Paper Source

The repository paper source has since been normalized to Type-A-style split TeX: `paper/main.tex`, `paper/abstract.tex`, `paper/macros.tex`, `paper/refs.bib`, `paper/sections/*.tex`, and `paper/appendices/*.tex`. The repository build now uses `pdflatex, bibtex, pdflatex, pdflatex` with `-jobname=Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer` and compiles cleanly to the titled PDF in `paper/`.
## Postscript: Visible PDF

Final paper red-team of the split source was completed on 2026-07-07. One S12 prose typo was fixed, the compressed certificate-dependent proof steps were expanded after external red-team review, the root split TeX files were synced into the repository package, and the visible 12-page reviewer PDF was generated at:

```text
paper/Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer.pdf
```

The final LaTeX log has no real LaTeX warnings, undefined references, overfull boxes, or fatal errors. BibTeX reports `warning$ -- 0`. The MiKTeX console still emits its environment-level update reminder, which is not a manuscript warning.

## Postscript: Red-Team Polish Pass

After external red-team review on 2026-07-07, the paper source was expanded at the compressed certificate-dependent steps: the H/K true-scalar links now point inline to their exact replay certificates, the ratio-induction paragraph identifies the recurrence/sign certificates, and the S12 tiler proof states the finite-window and symbolic tail checks more explicitly. The targeted verifiers S11A-H, S11A-transfer, S11A-K local/global, S8, S9C, S12A, and S12B were rerun and passed. `SHA256SUMS.txt` was regenerated with LF line endings.
## Boundary

This directory is the final local public-repo package. After human PDF review, it is ready to initialize, push, and tag as the first public release. The heavy rational audits remain an optional replay tier, not a blocker for the scoped theorem text.

## Next Task

Human review of paper/Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer.pdf; then apply any reader-facing edits and tag the first public release.
