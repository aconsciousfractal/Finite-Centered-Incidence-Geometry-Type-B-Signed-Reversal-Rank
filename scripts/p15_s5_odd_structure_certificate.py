# -*- coding: utf-8 -*-
"""
P15-S5 odd-n structure and reduction certificate builder.

This gate does not close the full odd-n theorem.  It records the exact odd
structure now available:

* V_ref for odd n has the split form
      M_ref = 2a P_plus_pairs + a_c E_center.
* The center scalar satisfies
      a_c(2m+1,k) = d_{k-1} - d_k,
  where d_p = |B~I_{2m}(p,p)|, and is positive by the even diagonal
  monotonicity consequence of the permanent-domination lemma.
* The non-center scalar a is positive for odd n>=5 and k>=2 by the trace
  reduction; the base n=3,k=1 is handled directly by exact rank/nonzero scalar.
* The remaining odd k=1, n>=5 V_ref residual is isolated as
      a(2m+1,1) = a(2m,1) + H_m.
  Positivity of H_m is replayed for finite range, but the recurrence certificate
  proving H_m>0 for all m is deferred to P15-S8.
* The odd pair-standard channel reduces exactly to lam_prime and mu; finite
  positivity is replayed for n=3,5,7, while the general proof remains open.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import factorial
from pathlib import Path

import p15_s3_exact_engine as s3

SCHEMA = "p15.s5.odd_structure_certificate.v1"


def frac_str(value: int | Fraction) -> str:
    value = value if isinstance(value, Fraction) else Fraction(value, 1)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def d_even(n_even: int, p: int, h_table: dict[tuple[int, int], int] | None = None) -> int:
    if p < 0:
        return 0
    if h_table is None:
        h_table = s3.h_even_table(n_even)
    return s3.size_even(n_even, p, p, h_table)


def odd_channel_exact_replay() -> dict:
    rows = []
    ok = True
    for n in (3, 5, 7):
        layers = s3.build_layers(n)
        for k in range(1, (n + 1) // 2 + 1):
            if (k, k) not in layers or not layers[(k, k)]:
                continue
            row = s3.channel_case(n, k, layers)
            a_value = row["a_actual"]
            ac_value = row["a_c_actual"]
            case_kind = "n3_base" if n == 3 else ("odd_k1_residual" if k == 1 else "odd_k_ge_2")
            a_condition_ok = (a_value != 0) if n == 3 else (a_value > 0)
            row_ok = (
                row["rank_ok"]
                and row["collapse_ok"]
                and row["a_c_closed_form_ok"]
                and row["odd_trace_ok"]
                and ac_value > 0
                and a_condition_ok
            )
            ok = ok and row_ok
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "case_kind": case_kind,
                    "size": row["size"],
                    "rank_M_ref": row["rank_M_ref"],
                    "rank_M_pair_std": row["rank_M_pair_std"],
                    "expected_rank_M_ref": row["expected_rank_M_ref"],
                    "expected_rank_M_pair_std": row["expected_rank_M_pair_std"],
                    "collapse_ok": row["collapse_ok"],
                    "a": a_value,
                    "a_nonzero": a_value != 0,
                    "a_positive": a_value > 0,
                    "a_c": ac_value,
                    "a_c_closed_form": row["a_c_closed_form"],
                    "a_c_positive": ac_value > 0,
                    "a_c_closed_form_ok": row["a_c_closed_form_ok"],
                    "odd_trace_ok": row["odd_trace_ok"],
                    "case_ok": row_ok,
                    "note": "n=3 has a=-1 but rank is correct because a is nonzero" if n == 3 else "finite replay supports the odd reduction",
                }
            )
    return {
        "id": "odd_channel_exact_replay",
        "status": "PASS" if ok else "FAIL",
        "range": "exact full enumeration for odd n=3,5,7 diagonal nonempty k>=1",
        "rows": rows,
    }


def odd_ac_closed_form_sweep(limit_m: int = 10) -> dict:
    rows = []
    ok = True
    empty_boundary_rows = []
    for m in range(1, limit_m + 1):
        n_even = 2 * m
        n_odd = 2 * m + 1
        h = s3.h_even_table(n_even)
        values = []
        for k in range(1, m + 2):
            d_prev = d_even(n_even, k - 1, h)
            d_curr = d_even(n_even, k, h)
            ac = d_prev - d_curr
            nonempty_by_center_scalar = ac > 0
            boundary_empty_top = ac == 0
            nonnegative_ok = ac >= 0
            two_step_bound_ok = (k * k * d_curr <= d_prev) if nonempty_by_center_scalar else True
            strict_decrease_ok = d_curr < d_prev if nonempty_by_center_scalar else True
            ok = ok and nonnegative_ok and two_step_bound_ok and strict_decrease_ok
            if boundary_empty_top:
                empty_boundary_rows.append({"n_odd": n_odd, "k": k})
            values.append(
                {
                    "k": k,
                    "d_k_minus_1": d_prev,
                    "d_k": d_curr,
                    "a_c": ac,
                    "nonempty_by_center_scalar": nonempty_by_center_scalar,
                    "boundary_empty_top_candidate": boundary_empty_top,
                    "a_c_positive_on_nonempty_row": nonempty_by_center_scalar,
                    "strict_decrease_ok_on_nonempty_row": strict_decrease_ok,
                    "d_k_times_k_squared_le_d_prev_on_nonempty_row": two_step_bound_ok,
                }
            )
        rows.append({"m": m, "n_odd": n_odd, "n_even_source": n_even, "values": values})
    return {
        "id": "odd_center_scalar_sweep",
        "status": "PASS" if ok else "FAIL",
        "range": f"odd n=3..{2 * limit_m + 1}, k=1..m+1 via even d_p counts; a_c>0 exactly on nonempty center-scalar rows",
        "rows": rows,
        "empty_boundary_rows": empty_boundary_rows,
        "formula": "a_c(2m+1,k)=d_{k-1}-d_k, d_p=|B~I_{2m}(p,p)|; theorem rows require nonempty B~I_{2m+1}(k,k)",
    }


# Univariate polynomial helpers for H_m.

def poly_mul(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int)
    for i, ci in a.items():
        for j, cj in b.items():
            out[i + j] += ci * cj
    return dict(out)


def poly_add(*polys: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int)
    for poly in polys:
        for i, coeff in poly.items():
            out[i] += coeff
    return dict(out)


def poly_scale(poly: dict[int, int], scalar: int) -> dict[int, int]:
    return {i: coeff * scalar for i, coeff in poly.items()}


def poly_pow(poly: dict[int, int], power: int) -> dict[int, int]:
    out = {0: 1}
    for _ in range(power):
        out = poly_mul(out, poly)
    return out


def L_func(poly: dict[int, int], n: int) -> int:
    return sum(coeff * (2 ** (n - r)) * factorial(n - r) for r, coeff in poly.items() if 0 <= r <= n)


A_POLY = {0: 1, 1: -4, 2: 2}


def B_m_agent1(m: int) -> int:
    return L_func(poly_mul({0: 1, 1: -1}, poly_pow(A_POLY, m - 1)), 2 * m - 1)


def S_m_agent1(m: int) -> int:
    # R_m(z,s)=Az^{m-2}(1+(z-2)s)(1-s)(1-2s), then Phi(R_m(1)-dz R_m(1)).
    def z_mul(a: dict[tuple[int, int], int], b: dict[tuple[int, int], int]) -> dict[tuple[int, int], int]:
        out: dict[tuple[int, int], int] = defaultdict(int)
        for (zi, si), ci in a.items():
            for (zj, sj), cj in b.items():
                out[(zi + zj, si + sj)] += ci * cj
        return dict(out)

    def z_pow(poly: dict[tuple[int, int], int], power: int) -> dict[tuple[int, int], int]:
        out = {(0, 0): 1}
        for _ in range(power):
            out = z_mul(out, poly)
        return out

    az = {
        (0, 0): 1,
        (1, 1): 2,
        (0, 1): -6,
        (2, 2): 1,
        (1, 2): -4,
        (0, 2): 5,
    }
    tail = {(0, 0): 1, (1, 1): 1, (0, 1): -2}
    one_minus_s = {(0, 0): 1, (0, 1): -1}
    one_minus_2s = {(0, 0): 1, (0, 1): -2}
    z_poly = z_mul(z_mul(z_mul(z_pow(az, m - 2), tail), one_minus_s), one_minus_2s)

    at_z1: dict[int, int] = defaultdict(int)
    dz_at_z1: dict[int, int] = defaultdict(int)
    for (z_degree, s_degree), coeff in z_poly.items():
        at_z1[s_degree] += coeff
        dz_at_z1[s_degree] += coeff * z_degree
    univariate = poly_add(dict(at_z1), poly_scale(dict(dz_at_z1), -1))
    return L_func(univariate, 2 * m - 1)


def H_m_agent1(m: int) -> int:
    return (2 * m - 2) * S_m_agent1(m) - B_m_agent1(m)


def B_m_agent2(m: int) -> int:
    return L_func(poly_mul({0: 1, 1: -1}, poly_pow(A_POLY, m - 1)), 2 * m - 1)


def L_m_agent2(m: int) -> int:
    cubic = {3: 2 * m, 2: -(4 * m + 2), 1: 2 * m + 2, 0: -1}
    poly = poly_mul(poly_mul(poly_mul({0: 1, 1: -1}, {0: -1, 1: 2}), poly_pow(A_POLY, m - 3)), cubic)
    return L_func(poly, 2 * m - 1)


def H_m_agent2(m: int) -> int:
    cubic = {3: 2 * m, 2: -(4 * m + 2), 1: 2 * m + 2, 0: -1}
    p_m = poly_add(poly_scale(poly_mul({0: -1, 1: 2}, cubic), 2 * m - 2), poly_scale(poly_pow(A_POLY, 2), -1))
    poly = poly_mul(poly_mul({0: 1, 1: -1}, poly_pow(A_POLY, m - 3)), p_m)
    return L_func(poly, 2 * m - 1)


# Truncated polynomial permanent for the odd k=1 ground residual.

def uv_t_add(a: dict[tuple[int, int], dict[int, int]], b: dict[tuple[int, int], dict[int, int]]) -> dict[tuple[int, int], dict[int, int]]:
    out: dict[tuple[int, int], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for poly in (a, b):
        for key, t_poly in poly.items():
            for t_degree, coeff in t_poly.items():
                out[key][t_degree] += coeff
    return {key: dict(value) for key, value in out.items()}


def uv_t_mul(a: dict[tuple[int, int], dict[int, int]], b: dict[tuple[int, int], dict[int, int]]) -> dict[tuple[int, int], dict[int, int]]:
    out: dict[tuple[int, int], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for (au, av), a_t in a.items():
        for (bu, bv), b_t in b.items():
            du, dv = au + bu, av + bv
            if du > 1 or dv > 1:
                continue
            for ta, ca in a_t.items():
                for tb, cb in b_t.items():
                    out[(du, dv)][ta + tb] += ca * cb
    return {key: dict(value) for key, value in out.items()}


def odd_layer11_negfix_poly(n: int) -> dict[int, int]:
    rho = s3.reversal(n)
    center = (n - 1) // 2

    def cell(i: int, j: int) -> dict[tuple[int, int], dict[int, int]]:
        if i == center and j == center:
            return {(1, 1): {0: 1}, (0, 0): {1: 1}}
        if j == i:
            return {(1, 0): {0: 1}, (0, 0): {1: 1}}
        if j == rho[i]:
            return {(0, 1): {0: 1}, (0, 0): {0: 1}}
        return {(0, 0): {0: 2}}

    total: dict[tuple[int, int], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for r in range(1, n + 1):
        for subset in combinations(range(n), r):
            prod = {(0, 0): {0: 1}}
            for i in range(n):
                row_sum: dict[tuple[int, int], dict[int, int]] = {}
                for j in subset:
                    row_sum = uv_t_add(row_sum, cell(i, j))
                prod = uv_t_mul(prod, row_sum)
                if not prod:
                    break
            sign = -1 if (n - r) % 2 else 1
            for key, t_poly in prod.items():
                for t_degree, coeff in t_poly.items():
                    total[key][t_degree] += sign * coeff
    return dict(total.get((1, 1), {}))


def even_a_k1(n_even: int, h_table: dict[tuple[int, int], int] | None = None) -> int:
    if h_table is None:
        h_table = s3.h_even_table(n_even)
    return s3.a_even_closed_form(n_even, 1, h_table)


def odd_k1_ground_residual(m: int) -> dict:
    n = 2 * m + 1
    h = s3.h_even_table(2 * m)
    q = odd_layer11_negfix_poly(n)
    size = sum(q.values())
    sum_neg_fixed = sum(t_degree * coeff for t_degree, coeff in q.items())
    ac = d_even(2 * m, 0, h) - d_even(2 * m, 1, h)
    na_odd = size - sum_neg_fixed - ac
    a_even = even_a_k1(2 * m, h)
    residual_c = na_odd - 2 * m * a_even
    return {
        "m": m,
        "n": n,
        "size": size,
        "sum_neg_fixed": sum_neg_fixed,
        "a_c": ac,
        "two_m_times_a_odd": na_odd,
        "a_odd": frac_str(Fraction(na_odd, 2 * m)),
        "a_even_2m_1": a_even,
        "residual_C": residual_c,
    }


def odd_k1_H_reduction(max_ground_m: int = 7, max_positive_m: int = 16) -> dict:
    ground_rows = []
    ok = True
    for m in range(2, max_ground_m + 1):
        ground = odd_k1_ground_residual(m)
        h1 = H_m_agent1(m)
        c1 = 2 * m * h1
        h2 = H_m_agent2(m) if m >= 3 else None
        c2 = 2 * m * h2 if h2 is not None else None
        match1 = c1 == ground["residual_C"]
        match2 = True if c2 is None else c2 == ground["residual_C"]
        identity_ok = ground["two_m_times_a_odd"] == 2 * m * (ground["a_even_2m_1"] + h1)
        ok = ok and match1 and match2 and identity_ok
        row = dict(ground)
        row.update(
            {
                "H_m_agent1": h1,
                "two_m_H_m_agent1": c1,
                "agent1_matches_ground_C": match1,
                "H_m_agent2": h2,
                "two_m_H_m_agent2": c2,
                "agent2_matches_ground_C": match2,
                "a_odd_equals_a_even_plus_H_m": identity_ok,
            }
        )
        ground_rows.append(row)

    positivity_rows = []
    all_positive = True
    for m in range(2, max_positive_m + 1):
        h = H_m_agent1(m)
        b = B_m_agent1(m)
        s_value = S_m_agent1(m)
        positive = h > 0
        all_positive = all_positive and positive
        positivity_rows.append(
            {
                "m": m,
                "H_m": h,
                "H_m_positive": positive,
                "B_m": b,
                "S_m": s_value,
                "ratio_string": f"{(2 * m - 2) * s_value}/{b}",
            }
        )
    ok = ok and all_positive
    return {
        "id": "odd_k1_H_reduction",
        "status": "PASS" if ok else "FAIL",
        "ground_range": f"m=2..{max_ground_m} exact truncated-permanent residual C_m",
        "positivity_range": f"m=2..{max_positive_m} H_m finite positivity replay",
        "identity": "a(2m+1,1)=a(2m,1)+H_m for m>=2; H_m>0 is the S8 residual positivity target",
        "ground_rows": ground_rows,
        "positivity_rows": positivity_rows,
        "open_dependency": "The all-m proof of H_m>0 still requires the P15-S8 recurrence/creative-telescoping certificate.",
    }


def pair_basis_matrix(matrix: list[list[int]], n: int) -> tuple[list[list[int]], list[tuple[int, int]], int]:
    rho = s3.reversal(n)
    center = (n - 1) // 2
    pairs = []
    seen: set[int] = set()
    for i in range(n):
        if i == center or i in seen:
            continue
        pairs.append((i, rho[i]))
        seen.add(i)
        seen.add(rho[i])

    basis_matrix: list[list[int]] = []
    for row_rep, _ in pairs:
        row = []
        for i, ri in pairs:
            row.append(matrix[row_rep][i] + matrix[row_rep][ri])
        row.append(matrix[row_rep][center])
        basis_matrix.append(row)
    center_row = []
    for i, ri in pairs:
        center_row.append(matrix[center][i] + matrix[center][ri])
    center_row.append(matrix[center][center])
    basis_matrix.append(center_row)
    return basis_matrix, pairs, center


def mat_vec(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def vec_equals_scaled(actual: list[int], scalar: Fraction, vector: list[int]) -> bool:
    return all(Fraction(value, 1) == scalar * vector[i] for i, value in enumerate(actual))


def odd_pair_standard_reduction() -> dict:
    rows = []
    ok = True
    for n in (3, 5, 7):
        layers = s3.build_layers(n)
        h_even = s3.h_even_table(n - 1)
        for k in range(1, (n + 1) // 2 + 1):
            x = layers.get((k, k), [])
            if not x:
                continue
            _, m_pair = s3.channel_matrices(x, n)
            basis_matrix, pairs, center = pair_basis_matrix(m_pair, n)
            m = (n - 1) // 2
            size = len(x)
            g_center = m_pair[center][center]
            g_formula = d_even(n - 1, k - 1, h_even) + d_even(n - 1, k, h_even)
            mu = Fraction(n * g_center - size, n - 1)
            mu_vector = [1 for _ in range(m)] + [-(n - 1)]
            mu_ok = vec_equals_scaled(mat_vec(basis_matrix, mu_vector), mu, mu_vector)

            lam_value = None
            lam_ok = True
            lam_positive = None
            if m >= 2:
                i0, ri0 = pairs[0]
                generic = [j for j in range(n) if j not in (i0, ri0, center)]
                lam_value = m_pair[i0][i0] + m_pair[ri0][i0] - 2 * m_pair[generic[0]][i0]
                for p in range(m - 1):
                    vector = [0 for _ in range(m + 1)]
                    vector[p] = 1
                    vector[m - 1] = -1
                    lam_ok = lam_ok and vec_equals_scaled(mat_vec(basis_matrix, vector), Fraction(lam_value, 1), vector)
                lam_positive = lam_value > 0

            ones_vector = [1 for _ in range(m)] + [1]
            ones_ok = vec_equals_scaled(mat_vec(basis_matrix, ones_vector), Fraction(size, 1), ones_vector)
            rank_pair_std = s3.rank_exact(s3.standard_projected_matrix(m_pair))
            expected_rank = m
            case_ok = (
                ones_ok
                and mu_ok
                and lam_ok
                and g_center == g_formula
                and mu > 0
                and (lam_positive is not False)
                and rank_pair_std == expected_rank
            )
            ok = ok and case_ok
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "size": size,
                    "rank_M_pair_std": rank_pair_std,
                    "expected_rank_M_pair_std": expected_rank,
                    "W_plus_basis_matrix": basis_matrix,
                    "g_center": g_center,
                    "g_center_formula_d_prev_plus_d_curr": g_formula,
                    "g_center_formula_ok": g_center == g_formula,
                    "lam_prime": lam_value,
                    "lam_prime_multiplicity": max(m - 1, 0),
                    "lam_prime_positive": lam_positive,
                    "lam_prime_eigenvector_check_ok": lam_ok,
                    "mu": frac_str(mu),
                    "mu_multiplicity": 1,
                    "mu_positive": mu > 0,
                    "mu_eigenvector_check_ok": mu_ok,
                    "all_ones_eigenvalue_ok": ones_ok,
                    "case_ok": case_ok,
                }
            )
    return {
        "id": "odd_pair_standard_reduction",
        "status": "PASS" if ok else "FAIL",
        "range": "exact full enumeration n=3,5,7; eigenvectors checked over Q/int arithmetic",
        "rows": rows,
        "open_dependency": "General positivity of lam_prime and mu remains a later proof obligation; S5 records the exact reduction.",
    }


def make_report(limit_m: int = 10, hm_ground_m: int = 7, hm_positive_m: int = 16) -> dict:
    checks = [
        odd_channel_exact_replay(),
        odd_ac_closed_form_sweep(limit_m),
        odd_k1_H_reduction(hm_ground_m, hm_positive_m),
        odd_pair_standard_reduction(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    proof_checklist = {
        "source_lock_S1": True,
        "prior_art_boundary_S2": True,
        "exact_engine_S3": True,
        "even_theorem_S4_available_for_even_a_and_d_counts": True,
        "odd_V_ref_split_form_recorded": True,
        "odd_n3_base_exact_nonzero_scalar": True,
        "odd_center_scalar_formula_and_positivity_recorded": True,
        "odd_noncenter_a_positive_for_k_ge_2_recorded": True,
        "odd_k1_reduced_to_H_m_residual_not_promoted": True,
        "odd_pair_standard_lam_prime_mu_reduction_recorded": True,
        "no_full_odd_theorem_promoted": True,
        "no_public_claim_promoted": True,
    }
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S5 odd-n structure",
        "status": "PASS" if not failed and all(proof_checklist.values()) else "FAIL",
        "decision": "S5 reduction gate passes; full odd-n theorem remains open pending H_m and pair-standard positivity proofs.",
        "failed_checks": failed,
        "scope": "odd n structure for the two ambient standard channels; no public theorem promotion",
        "checks": checks,
        "proof_checklist": proof_checklist,
        "closed_at_S5": [
            "Odd V_ref split form is recorded and replayed for n=3,5,7.",
            "Odd center scalar a_c=d_{k-1}-d_k is recorded and swept by exact even counts, with empty top rows fenced.",
            "Odd n=3,k=1 base is exact rank-correct even though a=-1.",
            "Odd k>=2 V_ref positivity proof is recorded as a reduction using permanent domination, center deletion, and add-center injection.",
            "Odd k=1 residual is isolated as a(2m+1,1)=a(2m,1)+H_m for m>=2.",
            "Odd pair-standard channel is reduced to exact lam_prime and mu eigenvalues.",
        ],
        "open_after_S5": [
            "P15-S8 must certify the H_m recurrence/creative-telescoping proof for all m.",
            "General positivity proofs for odd pair-standard lam_prime and mu remain open.",
            "No Type-D, C2xC2 square, full B_n fingerprint, or public theorem claim is promoted at S5.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S5 odd-n structure certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"], "-", check.get("range", check.get("ground_range", "")))
        if check["id"] == "odd_channel_exact_replay":
            for row in check["rows"]:
                print(
                    "  n={n} k={k} kind={case_kind} size={size} ranks=({rank_M_ref},{rank_M_pair_std}) "
                    "a={a} ac={a_c} collapse={collapse_ok} ok={case_ok}".format(**row)
                )
        elif check["id"] == "odd_center_scalar_sweep":
            for row in check["rows"]:
                vals = ", ".join("k=%d:ac=%d" % (v["k"], v["a_c"]) for v in row["values"])
                print("  n=%d %s" % (row["n_odd"], vals))
        elif check["id"] == "odd_k1_H_reduction":
            for row in check["ground_rows"]:
                print(
                    "  m={m} n={n} a_odd={a_odd} a_even={a_even_2m_1} H={H_m_agent1} "
                    "C={residual_C} matches={agent1_matches_ground_C}".format(**row)
                )
            last = check["positivity_rows"][-1]
            print("  H_m positive through m=%d; H_m=%d" % (last["m"], last["H_m"]))
        elif check["id"] == "odd_pair_standard_reduction":
            for row in check["rows"]:
                print(
                    "  n={n} k={k} rank={rank_M_pair_std} lam={lam_prime} mu={mu} "
                    "g_ok={g_center_formula_ok} eig_ok={lam_prime_eigenvector_check_ok}/{mu_eigenvector_check_ok}".format(**row)
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S5 odd-n structure certificate builder")
    parser.add_argument("--write-json", help="write deterministic JSON certificate to this path")
    parser.add_argument("--limit-m", type=int, default=10, help="maximum m for a_c sweep over odd n=2m+1")
    parser.add_argument("--hm-ground-m", type=int, default=7, help="maximum m for exact H_m residual ground replay")
    parser.add_argument("--hm-positive-m", type=int, default=16, help="maximum m for finite H_m positivity replay")
    args = parser.parse_args()
    report = make_report(args.limit_m, args.hm_ground_m, args.hm_positive_m)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
