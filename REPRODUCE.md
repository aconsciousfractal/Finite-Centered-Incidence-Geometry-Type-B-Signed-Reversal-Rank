# Reproduce

Run commands from the repository root. Use `python -B` so Python does not write `__pycache__` directories.

## Environment

```text
Python 3.11+
numpy
sympy
```

Install dependencies with:

```text
pip install -r requirements.txt
```

## Core Theorem Replay

```text
python -B scripts/p15_s4_even_theorem_certificate.py --write-json artifacts/certified/P15_S4_EVEN_N_THEOREM_CERTIFICATE.json
python -B scripts/p15_s5_odd_structure_certificate.py --write-json artifacts/certified/P15_S5_ODD_STRUCTURE_CERTIFICATE.json
python -B scripts/p15_s7_independent_red_team.py --write-json artifacts/certified/P15_S7_INDEPENDENT_RED_TEAM_CERTIFICATE.json
python -B scripts/p15_s8_hm_recurrence_certificate.py --write-json artifacts/certified/P15_S8_HM_RECURRENCE_CERTIFICATE.json
python -B scripts/p15_s9_odd_pair_reduction_certificate.py --write-json artifacts/certified/P15_S9_ODD_PAIR_REDUCTION_CERTIFICATE.json
python -B scripts/p15_s9c_km_recurrence_certificate.py --write-json artifacts/certified/P15_S9C_KM_RECURRENCE_CERTIFICATE.json
```

Expected result: every script exits with status `0` and reports `PASS` internally.

## Scalar-Link Replay

```text
python -B scripts/p15_s11a_h_true_scalar_edge_insertion_verifier.py --write-json artifacts/certified/P15_S11A_H_TRUE_SCALAR_EDGE_INSERTION_VERIFIER.json --write-md artifacts/certified/P15_S11A_H_TRUE_SCALAR_EDGE_INSERTION_VERIFIER.md
python -B scripts/p15_s11a_scalar_link_symbolic_verifier.py --write-json artifacts/certified/P15_S11A_SCALAR_LINK_SYMBOLIC_VERIFIER.json
python -B scripts/p15_s11a_k_true_scalar_bivariate_verifier.py --write-json artifacts/certified/P15_S11A_K_TRUE_SCALAR_BIVARIATE_VERIFIER.json
python -B scripts/p15_s11a_k_global_product_phi_certificate.py --write-json artifacts/certified/P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE.json --write-md artifacts/certified/P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE.md
python -B scripts/p15_s11_raw_rank_smoke_check.py --write-json artifacts/certified/P15_S11_RAW_RANK_SMOKE_CHECK_CERTIFICATE.json --write-md artifacts/certified/P15_S11_RAW_RANK_SMOKE_CHECK_CERTIFICATE.md
```

## Separator Add-On Replay

```text
python -B scripts/p15_s12_gate_a_signed_poset_certificate.py --write-json artifacts/certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.json --write-md artifacts/certified/P15_S12_GATE_A_SIGNED_POSET_CERTIFICATE.md
python -B scripts/p15_s12_gate_b_tiler_certificate.py --max-tail-probe-n 120 --write-json artifacts/certified/P15_S12_GATE_B_TILER_CERTIFICATE.json --write-md artifacts/certified/P15_S12_GATE_B_TILER_CERTIFICATE.md
```

## Type-D Sublayer Replay

These certify the type-D parity sublayer `D~I_n(k,k) = B~I_n(k,k) ∩ D_n`. The Type-D
scripts use the flag convention `--out <json>` / `--no-write` (not `--write-json`), report
a top-level `status` field, and exit `0` on PASS. They need only the Python standard
library; `p15d_g22` additionally imports the shared engine from `scripts/`.

```text
python -B type_d/p15d_g15_C_closure_verify.py --out type_d/P15D_G15_C_CLOSURE.json
python -B type_d/p15d_g20_oatail_symbolic_verify.py --out type_d/P15D_G20_OATAIL_SYMBOLIC.json
python -B type_d/p15d_g22_psdom_typeD_trace_verify.py --out type_d/P15D_G22_PSDOM_TYPED_TRACE.json
python -B type_d/p15d_g16_pairstd_closure_verify.py --out type_d/P15D_G16_PAIRSTD_CLOSURE.json
```

`p15d_g16_pairstd_closure_verify.py` recomputes exact permanents up to `n = 13` and takes
several minutes; the other three finish in well under a minute. See `type_d/README.md`.

## Optional Full-Fingerprint Appendix Replay

This tier is finite `n=8` appendix evidence only. It is not required for the main two-standard-channel theorem.

```text
python -B scripts/p15_full_bn_zero_family_group_algebra_certificate.py --write-json artifacts/certified/P15_FULL_BN_ZERO_FAMILY_GROUP_ALGEBRA_CERTIFICATE.json --write-md artifacts/certified/P15_FULL_BN_ZERO_FAMILY_GROUP_ALGEBRA_CERTIFICATE.md
python -B scripts/p15_full_bn_pure_plus_johnson_degree_drop.py --write-json artifacts/certified/P15_FULL_BN_PURE_PLUS_JOHNSON_DEGREE_DROP_CERTIFICATE.json --write-md artifacts/certified/P15_FULL_BN_PURE_PLUS_JOHNSON_DEGREE_DROP_CERTIFICATE.md
python -B scripts/p15_full_bn_pure_plus_zero_certificate.py --write-json artifacts/certified/P15_FULL_BN_PURE_PLUS_ZERO_CERTIFICATE.json --write-md artifacts/certified/P15_FULL_BN_PURE_PLUS_ZERO_CERTIFICATE.md
python -B scripts/p15_full_bn_lambda41_k3_plus_standard_certificate.py --write-json artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.json --write-md artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.md
```

The exact rational `T_1/T_2/T_4` audit is heavier. Run it family by family:

```text
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T1 --write-json artifacts/certified/P15_FULL_BN_T1_RATIONAL_DIRECT_RANK_AUDIT.json --write-md artifacts/certified/P15_FULL_BN_T1_RATIONAL_DIRECT_RANK_AUDIT.md
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T2 --write-json artifacts/certified/P15_FULL_BN_T2_RATIONAL_DIRECT_RANK_AUDIT.json --write-md artifacts/certified/P15_FULL_BN_T2_RATIONAL_DIRECT_RANK_AUDIT.md
python -B scripts/p15_full_bn_t124_rational_direct_rank_audit.py --family T4 --write-json artifacts/certified/P15_FULL_BN_T4_RATIONAL_DIRECT_RANK_AUDIT.json --write-md artifacts/certified/P15_FULL_BN_T4_RATIONAL_DIRECT_RANK_AUDIT.md
```


## Optional Citation Metadata Check

`CITATION.cff` is valid YAML and includes the repository citation plus the preferred article citation. If `cffconvert` is installed, run:

```text
cffconvert --validate -i CITATION.cff
```

The code license is MIT; paper/prose/PDF artifacts are CC-BY-4.0 as specified in `LICENSE`.

## Paper Build

From the repository root:

```text
cd paper
pdflatex -interaction=nonstopmode -halt-on-error -jobname=Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer main.tex
bibtex Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer
pdflatex -interaction=nonstopmode -halt-on-error -jobname=Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer main.tex
pdflatex -interaction=nonstopmode -halt-on-error -jobname=Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer main.tex
```

Expected PDF: `paper/Exact_Standard-Channel_Rank_of_a_Signed-Reversal_Incidence_Layer.pdf`. Generated LaTeX byproducts are ignored by Git by default; the titled PDF is the reviewer preview artifact.
