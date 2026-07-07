# -*- coding: utf-8 -*-
"""
P15-S9 odd pair-standard reduction certificate.

This gate closes the elementary odd pair-standard positivity reductions after S8:

* mu is positive for all odd nonempty diagonal rows because
      mu(2m+1,k) = (d_{k-1}-d_k) + 4*a_even(2m,k).
* lambda_prime is positive for all odd nonempty diagonal rows with k>=2 by
  the trace identity
      (m-1)*lambda_prime = (k-1)*|X| + NegFixSum - mu,
  the bound mu<=|X|, and an explicit k=2 negative-fixed witness.
* The k=1 odd pair-standard residue is isolated as
      lambda_prime(2m+1,1)>0, equivalently K_m=(m-1)*lambda_prime>0.
  It is closed downstream by the P15-S9C K_m certificate.

No recurrence search is performed here. The k=1 formula is recorded as the
S9C certificate target; current governance uses S9C for the all-m proof.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from math import factorial
from pathlib import Path

import p15_s3_exact_engine as s3

SCHEMA = "p15.s9.odd_pair_reduction_certificate.v1"


def frac_str(value: int | Fraction) -> str:
    value = value if isinstance(value, Fraction) else Fraction(value, 1)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


# ---------------------------------------------------------------------------
# Small exact polynomial algebra for the mu bridge identity.

Monomial = tuple[int, int, int, int, int]  # N, K, Dprev, D, E
Poly = dict[Monomial, int]


def p_clean(poly: Poly) -> Poly:
    return {key: value for key, value in poly.items() if value}


def p_const(value: int) -> Poly:
    return {} if value == 0 else {(0, 0, 0, 0, 0): value}


def p_var(index: int) -> Poly:
    exp = [0, 0, 0, 0, 0]
    exp[index] = 1
    return {tuple(exp): 1}  # type: ignore[return-value]


def p_add(*polys: Poly) -> Poly:
    out: Poly = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = out.get(key, 0) + value
    return p_clean(out)


def p_neg(poly: Poly) -> Poly:
    return {key: -value for key, value in poly.items()}


def p_sub(a: Poly, b: Poly) -> Poly:
    return p_add(a, p_neg(b))


def p_scale(poly: Poly, scalar: int) -> Poly:
    return p_clean({key: scalar * value for key, value in poly.items()})


def p_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for exp_a, coeff_a in a.items():
        for exp_b, coeff_b in b.items():
            key = tuple(exp_a[i] + exp_b[i] for i in range(5))
            out[key] = out.get(key, 0) + coeff_a * coeff_b
    return p_clean(out)


def p_format(poly: Poly) -> list[dict[str, int | list[int]]]:
    return [
        {"exponents_N_K_Dprev_D_E": list(key), "coeff": value}
        for key, value in sorted(poly.items())
    ]


N_VAR = p_var(0)
K_VAR = p_var(1)
DPREV_VAR = p_var(2)
D_VAR = p_var(3)
E_VAR = p_var(4)


def mu_bridge_symbolic_identity() -> dict:
    n_plus_1 = p_add(N_VAR, p_const(1))
    dsum = p_add(DPREV_VAR, D_VAR)
    center_nonfixed = p_scale(
        p_add(
            p_mul(p_sub(N_VAR, p_scale(K_VAR, 2)), D_VAR),
            p_scale(p_mul(p_add(K_VAR, p_const(1)), E_VAR), 2),
        ),
        2,
    )
    size_expr = p_add(dsum, center_nonfixed)
    n_mu_from_definition = p_sub(p_mul(n_plus_1, dsum), size_expr)
    n_mu_bridge = p_add(
        p_mul(N_VAR, p_sub(DPREV_VAR, D_VAR)),
        p_scale(p_sub(p_mul(K_VAR, D_VAR), p_mul(p_add(K_VAR, p_const(1)), E_VAR)), 4),
    )
    residual = p_sub(n_mu_from_definition, n_mu_bridge)
    return {
        "id": "mu_bridge_symbolic_identity",
        "status": "PASS" if not residual else "FAIL",
        "identity": "mu=(d_{k-1}-d_k)+(4/N)*(k*d_k-(k+1)*e_k)=a_c+4*a_even",
        "size_formula": "|X|=d_{k-1}+d_k+2*((N-2k)*d_k+2*(k+1)*e_k)",
        "variables": ["N", "k", "d_{k-1}", "d_k", "e_k"],
        "residual_terms": p_format(residual),
    }


# ---------------------------------------------------------------------------
# Closed-form sweeps for mu.


def even_diag(n_even: int, k: int, h_table: dict[tuple[int, int], int]) -> int:
    if k < 0:
        return 0
    return s3.size_even(n_even, k, k, h_table)


def mu_bridge_closed_form_sweep(max_m: int) -> dict:
    rows = []
    ok = True
    for m in range(1, max_m + 1):
        n_odd = 2 * m + 1
        n_even = 2 * m
        h = s3.h_even_table(n_even)
        values = []
        for k in range(1, m + 2):
            d_prev = even_diag(n_even, k - 1, h)
            d_curr = even_diag(n_even, k, h)
            e_k = s3.size_even(n_even, k + 1, k, h)
            size_formula = d_prev + d_curr + 2 * ((n_even - 2 * k) * d_curr + 2 * (k + 1) * e_k)
            a_even = Fraction(k * d_curr - (k + 1) * e_k, n_even)
            mu_from_definition = Fraction((n_even + 1) * (d_prev + d_curr) - size_formula, n_even)
            mu_bridge = Fraction(d_prev - d_curr) + 4 * a_even
            nonempty = size_formula > 0
            mu_le_size = mu_from_definition <= size_formula if nonempty else True
            row_ok = mu_from_definition == mu_bridge and (mu_from_definition > 0 if nonempty else True) and mu_le_size
            ok = ok and row_ok
            values.append(
                {
                    "k": k,
                    "nonempty_by_size_formula": nonempty,
                    "size_formula": size_formula,
                    "d_k_minus_1": d_prev,
                    "d_k": d_curr,
                    "e_k": e_k,
                    "a_c": d_prev - d_curr,
                    "a_even": frac_str(a_even),
                    "mu_from_definition": frac_str(mu_from_definition),
                    "mu_bridge": frac_str(mu_bridge),
                    "mu_positive_on_nonempty": (mu_from_definition > 0) if nonempty else None,
                    "mu_le_size_on_nonempty": mu_le_size if nonempty else None,
                    "case_ok": row_ok,
                }
            )
        rows.append({"m": m, "n_odd": n_odd, "n_even_source": n_even, "values": values})
    return {
        "id": "mu_bridge_closed_form_sweep",
        "status": "PASS" if ok else "FAIL",
        "range": f"odd n=3..{2 * max_m + 1} using even closed counts only",
        "proof_dependencies": [
            "S5: a_c(2m+1,k)=d_{k-1}-d_k is positive on nonempty odd rows",
            "S4: a_even(2m,k)>0 on nonempty even rows; top odd rows are covered by a_c",
            "S3/S4 trace identity: 2m*a_even(2m,k)=k*d_k-(k+1)*e_k",
        ],
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# Exact finite replay for the lambda trace identity.


def odd_pairs_and_center(n: int) -> tuple[list[tuple[int, int]], int]:
    rho = s3.reversal(n)
    center = (n - 1) // 2
    seen: set[int] = {center}
    pairs = []
    for i in range(n):
        if i in seen:
            continue
        pairs.append((i, rho[i]))
        seen.add(i)
        seen.add(rho[i])
    return pairs, center


def lambda_trace_finite_replay() -> dict:
    rows = []
    ok = True
    for n in (3, 5, 7):
        layers = s3.build_layers(n)
        pairs, center = odd_pairs_and_center(n)
        m = (n - 1) // 2
        for k in range(1, m + 2):
            layer = layers.get((k, k), [])
            if not layer:
                continue
            _, m_pair = s3.channel_matrices(layer, n)
            size = len(layer)
            g_center = m_pair[center][center]
            mu = Fraction(n * g_center - size, n - 1)
            trace_pair = sum(m_pair[i][i] for i in range(n))
            fixed_sum = sum(s3.fixed_count(pi) for pi, _eps in layer)
            neg_fixed_sum = sum(s3.neg_fixed_count(pi, eps) for pi, eps in layer)
            positive_fixed_identity = k * size + neg_fixed_sum
            lambda_prime: int | None = None
            lambda_trace_lhs = Fraction(0, 1)
            lambda_trace_rhs = Fraction((k - 1) * size + neg_fixed_sum, 1) - mu
            if m >= 2:
                i0, ri0 = pairs[0]
                generic = [j for j in range(n) if j not in (i0, ri0, center)]
                lambda_prime = m_pair[i0][i0] + m_pair[ri0][i0] - 2 * m_pair[generic[0]][i0]
                lambda_trace_lhs = Fraction(m - 1, 1) * lambda_prime
            trace_decomposition = Fraction(size, 1) + lambda_trace_lhs + mu
            fixed_trace_ok = trace_pair == fixed_sum == positive_fixed_identity
            lambda_formula_ok = lambda_trace_lhs == lambda_trace_rhs
            trace_decomposition_ok = Fraction(trace_pair, 1) == trace_decomposition
            mu_le_size = mu <= size
            row_ok = fixed_trace_ok and lambda_formula_ok and trace_decomposition_ok and mu_le_size
            ok = ok and row_ok
            rows.append(
                {
                    "n": n,
                    "m": m,
                    "k": k,
                    "size": size,
                    "g_center": g_center,
                    "mu": frac_str(mu),
                    "mu_le_size": mu_le_size,
                    "lambda_prime": lambda_prime,
                    "neg_fixed_sum": neg_fixed_sum,
                    "trace_pair": trace_pair,
                    "fixed_sum": fixed_sum,
                    "k_size_plus_neg_fixed_sum": positive_fixed_identity,
                    "fixed_trace_ok": fixed_trace_ok,
                    "trace_decomposition": frac_str(trace_decomposition),
                    "trace_decomposition_ok": trace_decomposition_ok,
                    "lambda_trace_lhs": frac_str(lambda_trace_lhs),
                    "lambda_trace_rhs": frac_str(lambda_trace_rhs),
                    "lambda_trace_identity_ok": lambda_formula_ok,
                    "case_ok": row_ok,
                }
            )
    return {
        "id": "lambda_trace_finite_replay",
        "status": "PASS" if ok else "FAIL",
        "range": "exact full enumeration for odd n=3,5,7 diagonal nonempty k>=1",
        "identity": "(m-1)*lambda_prime=(k-1)*|X|+NegFixSum-mu",
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# k>=2 closure logic and constructive k=2 witness.


def is_perm(pi: tuple[int, ...]) -> bool:
    return sorted(pi) == list(range(len(pi)))


def k2_negative_fixed_witness(n: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if n < 5 or n % 2 == 0:
        raise ValueError("k=2 witness is for odd n>=5")
    c = (n - 1) // 2
    pi = list(range(n))
    eps = [-1 for _ in range(n)]
    eps[c] = 1
    eps[0] = 1
    pi[1] = n - 2
    pi[n - 2] = 1
    eps[1] = 1
    eps[n - 2] = -1
    return tuple(pi), tuple(eps)


def lambda_k_ge_2_closure(max_m: int) -> dict:
    witness_rows = []
    ok = True
    for m in range(2, max_m + 1):
        n = 2 * m + 1
        pi, eps = k2_negative_fixed_witness(n)
        hits = s3.positive_hit_counts(pi, eps, s3.reversal(n))
        neg_fixed = s3.neg_fixed_count(pi, eps)
        row_ok = is_perm(pi) and hits == (2, 2) and neg_fixed > 0
        ok = ok and row_ok
        witness_rows.append(
            {
                "n": n,
                "m": m,
                "positive_hit_counts": list(hits),
                "neg_fixed_count": neg_fixed,
                "permutation_valid": is_perm(pi),
                "case_ok": row_ok,
            }
        )
    return {
        "id": "lambda_k_ge_2_closure",
        "status": "PASS" if ok else "FAIL",
        "logical_closure": [
            "Trace identity gives (m-1)*lambda_prime=(k-1)*|X|+NegFixSum-mu.",
            "The center count g_center<=|X| gives mu=((2m+1)*g_center-|X|)/(2m)<=|X|.",
            "For k>=3, the right side is at least (k-2)*|X|+NegFixSum>0 on nonempty rows.",
            "For k=2, the explicit witness below proves NegFixSum>0 for every odd n>=5; hence |X|+NegFixSum-mu>=NegFixSum>0.",
            "The m=1 case n=3 has no lambda_prime eigenspace, so no k>=2 lambda_prime claim is needed there.",
        ],
        "witness_range": f"odd n=5..{2 * max_m + 1}",
        "witness_pattern": "center and 0 are positive fixed; 1 maps positively to rho(1); rho(1) maps back negatively; all remaining default fixed entries are negative",
        "witness_rows": witness_rows,
    }


# ---------------------------------------------------------------------------
# k=1 residue formula target.

TriKey = tuple[int, int, int]  # u degree, v degree, t degree
TriPoly = dict[TriKey, int]
UniPoly = dict[int, int]


ONE_TRI: TriPoly = {(0, 0, 0): 1}
U_TRI: TriPoly = {(1, 0, 0): 1}
V_TRI: TriPoly = {(0, 1, 0): 1}
T_TRI: TriPoly = {(0, 0, 1): 1}


def tri_add(*polys: TriPoly) -> TriPoly:
    out: dict[TriKey, int] = defaultdict(int)
    for poly in polys:
        for key, coeff in poly.items():
            out[key] += coeff
    return {key: coeff for key, coeff in out.items() if coeff}


def tri_neg(poly: TriPoly) -> TriPoly:
    return {key: -coeff for key, coeff in poly.items()}


def tri_sub(a: TriPoly, b: TriPoly) -> TriPoly:
    return tri_add(a, tri_neg(b))


def tri_scale(poly: TriPoly, scalar: int) -> TriPoly:
    return {key: scalar * coeff for key, coeff in poly.items() if scalar * coeff}


def tri_mul(a: TriPoly, b: TriPoly) -> TriPoly:
    out: dict[TriKey, int] = defaultdict(int)
    for (au, av, at), ac in a.items():
        for (bu, bv, bt), bc in b.items():
            out[(au + bu, av + bv, at + bt)] += ac * bc
    return {key: coeff for key, coeff in out.items() if coeff}


def tri_pow(poly: TriPoly, exponent: int) -> TriPoly:
    out = ONE_TRI
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = tri_mul(out, base)
        base = tri_mul(base, base)
        n >>= 1
    return out


def phi_factorial(poly: UniPoly, n: int) -> int:
    total = 0
    for degree, coeff in poly.items():
        if 0 <= degree <= n:
            total += coeff * (2 ** (n - degree)) * factorial(n - degree)
    return total


def lambda_k1_bivariate_formula(m: int) -> tuple[int, UniPoly]:
    if m < 2:
        raise ValueError("lambda k=1 formula is used for m>=2")
    x = tri_sub(U_TRI, ONE_TRI)
    y = tri_sub(V_TRI, ONE_TRI)
    z = tri_sub(tri_mul(U_TRI, V_TRI), ONE_TRI)
    a_poly = tri_add(
        ONE_TRI,
        tri_scale(tri_mul(tri_add(x, y), T_TRI), 2),
        tri_mul(tri_add(tri_mul(x, x), tri_mul(y, y)), tri_mul(T_TRI, T_TRI)),
    )
    c_poly = tri_add(ONE_TRI, tri_mul(z, T_TRI))
    sx_poly = tri_add(ONE_TRI, tri_mul(x, T_TRI))
    sy_poly = tri_add(ONE_TRI, tri_mul(y, T_TRI))
    h_poly = tri_add(ONE_TRI, tri_mul(tri_add(x, y), T_TRI))
    bracket = tri_sub(
        tri_mul(
            tri_pow(a_poly, m - 1),
            tri_add(tri_mul(tri_add(ONE_TRI, U_TRI), sx_poly), tri_mul(tri_add(ONE_TRI, V_TRI), sy_poly)),
        ),
        tri_scale(tri_mul(tri_pow(a_poly, m - 2), tri_mul(h_poly, h_poly)), 4),
    )
    expr = tri_mul(c_poly, bracket)
    coeff_uv: dict[int, int] = defaultdict(int)
    for (u_degree, v_degree, t_degree), coeff in expr.items():
        if u_degree == 1 and v_degree == 1:
            coeff_uv[t_degree] += coeff
    coeff_uv = {degree: coeff for degree, coeff in coeff_uv.items() if coeff}
    return phi_factorial(coeff_uv, 2 * m), coeff_uv


def uni_mul(a: UniPoly, b: UniPoly) -> UniPoly:
    out: dict[int, int] = defaultdict(int)
    for i, ai in a.items():
        for j, bj in b.items():
            out[i + j] += ai * bj
    return {degree: coeff for degree, coeff in out.items() if coeff}


def uni_pow(poly: UniPoly, exponent: int) -> UniPoly:
    out: UniPoly = {0: 1}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = uni_mul(out, base)
        base = uni_mul(base, base)
        n >>= 1
    return out


def lambda_k1_univariate_formula(m: int) -> tuple[int, UniPoly]:
    if m < 4:
        raise ValueError("univariate residue formula is stated for m>=4")
    coeff: dict[int, int] = defaultdict(int)

    def put(degree: int, value: int) -> None:
        coeff[degree] += value

    put(7, 8 * m * m - 24 * m + 8)
    put(6, -16 * m * m + 56 * m)
    put(5, -12 * m * m - 44 * m + 12)
    put(4, 56 * m * m + 24 * m - 92)
    put(3, -56 * m * m - 24 * m + 122)
    put(2, 24 * m * m + 18 * m - 68)
    put(1, -4 * m * m - 8 * m + 19)
    put(0, 2 * m - 3)
    l_poly = {degree + 1: 2 * value for degree, value in coeff.items() if value}
    b_poly = {0: 1, 1: -4, 2: 2}
    residue_poly = uni_mul(uni_pow(b_poly, m - 4), l_poly)
    return phi_factorial(residue_poly, 2 * m), residue_poly


def lambda_k1_residue_target(max_m: int) -> dict:
    rows = []
    ok = True
    for m in range(2, max_m + 1):
        lambda_bivar, coeff_uv = lambda_k1_bivariate_formula(m)
        lambda_univar = None
        univar_match = None
        univar_degree = None
        if m >= 4:
            lambda_univar, residue_poly = lambda_k1_univariate_formula(m)
            univar_match = lambda_univar == lambda_bivar
            univar_degree = max(residue_poly) if residue_poly else -1
            ok = ok and univar_match
        rows.append(
            {
                "m": m,
                "n": 2 * m + 1,
                "lambda_prime_k1_bivariate": lambda_bivar,
                "K_m": (m - 1) * lambda_bivar,
                "uv_t_degree": max(coeff_uv) if coeff_uv else -1,
                "uv_t_terms": len(coeff_uv),
                "lambda_prime_k1_univariate": lambda_univar,
                "univariate_degree": univar_degree,
                "univariate_matches_bivariate": univar_match,
            }
        )
    known = {2: 28, 3: 7320, 4: 1747216}
    known_ok = all(next(row for row in rows if row["m"] == m)["lambda_prime_k1_bivariate"] == value for m, value in known.items())
    ok = ok and known_ok
    return {
        "id": "lambda_k1_residue_target",
        "status": "PASS" if ok else "FAIL",
        "range": f"formula values m=2..{max_m}; univariate check for m>=4",
        "known_exact_values_checked": known,
        "known_exact_values_ok": known_ok,
        "formula_status": "recorded_as_S9C_certificate_target; downstream all-m positivity is certified by P15-S9C",
        "s9c_target": "K_m=(m-1)*lambda_prime(2m+1,1)>0 for all m>=2; certified downstream by P15-S9C.",
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# Report and CLI.


def make_report(max_m: int, max_formula_m: int) -> dict:
    checks = [
        mu_bridge_symbolic_identity(),
        mu_bridge_closed_form_sweep(max_m),
        lambda_trace_finite_replay(),
        lambda_k_ge_2_closure(max_m),
        lambda_k1_residue_target(max_formula_m),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    proof_checklist = {
        "source_lock_S1": True,
        "exact_engine_S3_available": True,
        "even_theorem_S4_available_for_a_even_positivity": True,
        "odd_structure_S5_available_for_center_scalar_and_pair_reduction": True,
        "H_m_S8_closed_odd_V_ref": True,
        "mu_bridge_identity_verified_symbolically": checks[0]["status"] == "PASS",
        "mu_positive_all_odd_nonempty_rows_by_S4_S5_bridge": checks[1]["status"] == "PASS",
        "lambda_trace_identity_replayed_exactly": checks[2]["status"] == "PASS",
        "lambda_prime_positive_for_k_ge_2_by_trace_and_witness": checks[3]["status"] == "PASS",
        "lambda_prime_k1_isolated_for_S9C": checks[4]["status"] == "PASS",
        "no_full_odd_pair_standard_theorem_promoted_by_S9_alone": True,
        "no_public_claim_promoted": True,
    }
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S9 odd pair-standard reduction",
        "status": "PASS" if not failed and all(proof_checklist.values()) else "FAIL",
        "decision": "S9A/S9B pass: mu>0 for all odd nonempty rows and lambda_prime>0 for k>=2. The k=1 residue is isolated as K_m and is closed downstream by the P15-S9C certificate.",
        "failed_checks": failed,
        "scope": "odd pair-standard positivity reduction for the two ambient standard channels; S9 isolates k=1, S9C supplies the downstream all-m proof",
        "checks": checks,
        "proof_checklist": proof_checklist,
        "closed_at_S9": [
            "Odd pair-standard mu positivity is reduced to a_c+4*a_even and closed from S4/S5.",
            "Odd pair-standard lambda_prime positivity is closed for every k>=2 nonempty row by trace, mu<=|X|, and a k=2 negative-fixed witness.",
            "The k=1 odd pair-standard residue is isolated as K_m=(m-1)*lambda_prime(2m+1,1)>0 for S9C.",
        ],
        "downstream_after_S9": [
            "P15-S9C certifies K_m>0 for all m>=2 by integration-by-parts recurrence and ratio induction.",
            "No Type-D, full B_n fingerprint, positive-only separator, or public theorem claim is promoted by this gate.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S9 odd pair-standard reduction certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        if check["id"] == "mu_bridge_symbolic_identity":
            print("  residual_terms:", len(check["residual_terms"]))
        elif check["id"] == "mu_bridge_closed_form_sweep":
            print("  range:", check["range"])
            last = check["rows"][-1]
            nonempty = [v for v in last["values"] if v["nonempty_by_size_formula"]]
            print("  last_n:", last["n_odd"], "nonempty_rows:", len(nonempty))
        elif check["id"] == "lambda_trace_finite_replay":
            print("  rows:", len(check["rows"]))
            for row in check["rows"]:
                if row["n"] in (5, 7) and row["k"] == 1:
                    print("  n={n} k={k} lambda={lambda_prime} mu={mu} trace_ok={case_ok}".format(**row))
        elif check["id"] == "lambda_k_ge_2_closure":
            print("  witness_range:", check["witness_range"])
        elif check["id"] == "lambda_k1_residue_target":
            for row in check["rows"][:4]:
                print("  m={m} lambda={lambda_prime_k1_bivariate} K={K_m}".format(**row))
            print("  s9c_target:", check["s9c_target"])


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S9 odd pair-standard reduction certificate")
    parser.add_argument("--max-m", type=int, default=10, help="closed-form sweep and witness range")
    parser.add_argument("--max-formula-m", type=int, default=10, help="k=1 residue formula value range")
    parser.add_argument("--write-json", help="write JSON certificate")
    args = parser.parse_args()
    report = make_report(args.max_m, args.max_formula_m)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())