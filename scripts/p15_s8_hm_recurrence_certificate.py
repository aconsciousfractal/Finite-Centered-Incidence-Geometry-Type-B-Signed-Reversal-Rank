"""
P15-S8 H_m all-m recurrence certificate.

This certificate is stdlib-only. It verifies a concrete integration-by-parts
certificate for the order-3 recurrence of H_m. No symbolic search is performed.

For m>=3 the factorial functional is written as an integral

    H_m = integral_0^infty exp(-x) * 2^m * B(x)^(m-3) * Cbar_m(x) dx,
    B(x)=2x^2-4x+1.

The verifier checks the polynomial identity

    B*S_x + ((m-3)*B_x - B)*S = B*R,

where R=sum_j p_j(m) 2^j B^j Cbar_{m+j}. This means

    exp(-x)2^m B^(m-3)R = d/dx(exp(-x)2^m B^(m-3)S),

and the boundary terms vanish because S(m,0)=0 and exp(-x) kills infinity.
Thus the recurrence holds for all m>=3. The m=2 edge case is checked exactly.
The already elementary ratio induction then gives H_m>0 for all m>=2.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

SCHEMA = "p15.s8.hm_recurrence_certificate.v1"

P0 = [0, 835520, 7491008, 20923648, 28054784, 20596224, 8488960, 1843200, 163840]
P1 = [222992, 2074768, 5792192, 7631808, 5444992, 2173440, 460800, 40960]
P2 = [102746, 853240, 1735544, 575920, -1983040, -2734976, -1518080, -399360, -40960]
P3 = [-693, -5400, -8908, 1016, 13600, 10880, 2560]
RECURRENCE = [P0, P1, P2, P3]
A_POLY = [1, -4, 2]

# S(x,m) coefficient rows. Outer index is x-degree 0..11; inner index is m-degree 0..9.
S_COEFFS = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-362992, -4040528, -14671936, -18324992, 3534272, 30320896, 31213568, 14950400, 3522560, 327680],
    [4865224, 53275200, 190348288, 238904704, -35037184, -374093056, -389178880, -187074560, -44113920, -4096000],
    [-26025328, -283732432, -1010455376, -1286180672, 125308864, 1908027264, 2014666240, 973854720, 229990400, 21299200],
    [71845984, 789779584, 2841212832, 3687982784, -222214272, -5293494272, -5674812928, -2759813120, -652451840, -60129280],
    [-112358848, -1260043584, -4632385408, -6141249920, 266773760, 8807709184, 9551533056, 4667371520, 1102233600, 100597760],
    [103573120, 1198066624, 4553341056, 6196859008, -186378752, -8970813952, -9845695488, -4832931840, -1136967680, -101908480],
    [-51236864, -634644928, -2591901504, -3729769472, -57376000, 5344483840, 6038028288, 2998190080, 702218240, 60948480],
    [5608064, 117790208, 682714240, 1210087936, 251643392, -1555555328, -1970872320, -1024983040, -241172480, -19660800],
    [7350784, 50694400, 37381888, -141334528, -213181440, 22554624, 214401024, 149749760, 38010880, 2621440],
    [-3636864, -30290688, -61600768, -16988160, 78794752, 94875648, 37109760, 1638400, -1310720, 0],
    [532224, 4502016, 9606144, 3780608, -10964992, -15319040, -7536640, -1310720, 0, 0],
]

Poly = dict[tuple[int, int], Fraction]


def clean(poly: Poly) -> Poly:
    return {key: value for key, value in poly.items() if value}


def const(value: int | Fraction) -> Poly:
    value = value if isinstance(value, Fraction) else Fraction(value, 1)
    return {} if value == 0 else {(0, 0): value}


def xpow(power: int, coeff: int | Fraction = 1) -> Poly:
    coeff = coeff if isinstance(coeff, Fraction) else Fraction(coeff, 1)
    return {} if coeff == 0 else {(power, 0): coeff}


def mpow(power: int, coeff: int | Fraction = 1) -> Poly:
    coeff = coeff if isinstance(coeff, Fraction) else Fraction(coeff, 1)
    return {} if coeff == 0 else {(0, power): coeff}


def add(*polys: Poly) -> Poly:
    out: Poly = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = out.get(key, Fraction(0, 1)) + value
    return clean(out)


def neg(poly: Poly) -> Poly:
    return {key: -value for key, value in poly.items()}


def sub(a: Poly, b: Poly) -> Poly:
    return add(a, neg(b))


def scale(poly: Poly, scalar: int | Fraction) -> Poly:
    scalar = scalar if isinstance(scalar, Fraction) else Fraction(scalar, 1)
    return clean({key: scalar * value for key, value in poly.items()})


def mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for (xi, mi), ci in a.items():
        for (xj, mj), cj in b.items():
            key = (xi + xj, mi + mj)
            out[key] = out.get(key, Fraction(0, 1)) + ci * cj
    return clean(out)


def pow_poly(poly: Poly, exponent: int) -> Poly:
    out = const(1)
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def deriv_x(poly: Poly) -> Poly:
    out: Poly = {}
    for (xdeg, mdeg), coeff in poly.items():
        if xdeg:
            out[(xdeg - 1, mdeg)] = out.get((xdeg - 1, mdeg), Fraction(0, 1)) + coeff * xdeg
    return clean(out)


def poly_m_from_coeffs(coeffs: list[int]) -> Poly:
    return clean({(0, degree): Fraction(coeff, 1) for degree, coeff in enumerate(coeffs) if coeff})


def poly_t_list_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ci in enumerate(a):
        if ci:
            for j, cj in enumerate(b):
                if cj:
                    out[i + j] += ci * cj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_t_list_pow(poly: list[int], exponent: int) -> list[int]:
    out = [1]
    for _ in range(exponent):
        out = poly_t_list_mul(out, poly)
    return out


def poly_t_list_add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i, coeff in enumerate(a):
        out[i] += coeff
    for i, coeff in enumerate(b):
        out[i] += coeff
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_t_list_scale(poly: list[int], scalar: int) -> list[int]:
    out = [scalar * coeff for coeff in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def hm_value(m: int) -> int:
    if m == 2:
        return 22
    if m < 2:
        raise ValueError("H_m is used only for m>=2")
    cubic = [-1, 2 * m + 2, -(4 * m + 2), 2 * m]
    term = poly_t_list_mul([-1, 2], cubic)
    p_m = poly_t_list_add(poly_t_list_scale(term, 2 * m - 2), poly_t_list_scale(poly_t_list_pow(A_POLY, 2), -1))
    f_poly = poly_t_list_mul(poly_t_list_mul([1, -1], poly_t_list_pow(A_POLY, m - 3)), p_m)
    n = 2 * m - 1
    return sum(coeff * (2 ** (n - r)) * factorial(n - r) for r, coeff in enumerate(f_poly) if 0 <= r <= n)


def polyval_univar(coeffs: list[int], m: int) -> int:
    total = 0
    power = 1
    for coeff in coeffs:
        total += coeff * power
        power *= m
    return total


def recurrence_residual(m: int) -> int:
    values = {m + j: hm_value(m + j) for j in range(4)}
    return sum(polyval_univar(poly, m) * values[m + j] for j, poly in enumerate(RECURRENCE))


def m_shift(shift: int) -> Poly:
    return add(mpow(1), const(shift))


def cbar(shift: int) -> Poly:
    """Cbar_{m+shift}(x)=C_{m+shift}(x)/2^(m+shift)."""
    mm = m_shift(shift)
    mm2 = mul(mm, mm)
    inner = const(0)
    terms = [
        (mul(mm2, xpow(3)), -8),
        (mul(mm2, xpow(2)), 16),
        (mul(mm2, xpow(1)), -10),
        (mm2, 2),
        (mul(mm, xpow(4)), 8),
        (mul(mm, xpow(3)), -8),
        (mul(mm, xpow(2)), -4),
        (mul(mm, xpow(1)), 6),
        (mm, -2),
        (xpow(4), -12),
        (xpow(3), 32),
        (xpow(2), -32),
        (xpow(1), 12),
        (const(1), -1),
    ]
    for poly, coeff in terms:
        inner = add(inner, scale(poly, coeff))
    return scale(mul(add(scale(xpow(1), 2), const(-1)), inner), Fraction(1, 2))


def s_certificate() -> Poly:
    out: Poly = {}
    for xdeg, row in enumerate(S_COEFFS):
        for mdeg, coeff in enumerate(row):
            if coeff:
                out[(xdeg, mdeg)] = Fraction(coeff, 1)
    return out


def b_poly() -> Poly:
    return add(scale(xpow(2), 2), scale(xpow(1), -4), const(1))


def recurrence_polynomial_r() -> Poly:
    b = b_poly()
    out = const(0)
    for shift, p_coeffs in enumerate(RECURRENCE):
        term = mul(poly_m_from_coeffs(p_coeffs), scale(pow_poly(b, shift), 2 ** shift))
        term = mul(term, cbar(shift))
        out = add(out, term)
    return out


def max_degrees(poly: Poly) -> dict[str, int]:
    if not poly:
        return {"x": -1, "m": -1, "terms": 0}
    return {
        "x": max(key[0] for key in poly),
        "m": max(key[1] for key in poly),
        "terms": len(poly),
    }


def compute_d_polynomial() -> list[int]:
    def trim_list(poly: list[int]) -> list[int]:
        out = list(poly)
        while len(out) > 1 and out[-1] == 0:
            out.pop()
        return out

    def add_list(a: list[int], b: list[int], scale_b: int = 1) -> list[int]:
        n = max(len(a), len(b))
        out = [0] * n
        for i, coeff in enumerate(a):
            out[i] += coeff
        for i, coeff in enumerate(b):
            out[i] += scale_b * coeff
        return trim_list(out)

    def mul_list(a: list[int], b: list[int]) -> list[int]:
        out = [0] * (len(a) + len(b) - 1)
        for i, ci in enumerate(a):
            for j, cj in enumerate(b):
                out[i + j] += ci * cj
        return trim_list(out)

    def scale_list(poly: list[int], scalar: int) -> list[int]:
        return trim_list([scalar * coeff for coeff in poly])

    m2 = [0, 0, 1]
    mp1_2 = [1, 2, 1]
    mp2_2 = [4, 4, 1]
    abs_p2 = scale_list(P2, -1)
    term1 = mul_list(mul_list(scale_list(m2, 16), mp1_2), abs_p2)
    term2 = mul_list(scale_list(m2, 4), P1)
    term4 = mul_list(mul_list(mul_list(scale_list(m2, 64), mp1_2), mp2_2), P3)
    return add_list(add_list(add_list(term1, term2, -1), P0, -1), term4, -1)


def dominance_certificate(coeffs: list[int], m_min: int, name: str) -> dict:
    negative_degrees = [i for i, c in enumerate(coeffs) if c < 0]
    candidates = []
    for d in range(1, len(coeffs)):
        if coeffs[d] <= 0:
            continue
        if negative_degrees and max(negative_degrees) >= d:
            continue
        if all(c >= 0 for c in coeffs[d:]):
            candidates.append(d)
    if not candidates:
        finite_ok = all(polyval_univar(coeffs, m) > 0 for m in range(m_min, 200))
        return {"name": name, "status": "FAIL", "finite_sanity_ok": finite_ok, "method": "no coefficient-domination split found"}
    d = candidates[0]
    neg_sum = sum(-c for c in coeffs[:d] if c < 0)
    c_d = coeffs[d]
    threshold = max(m_min, neg_sum // c_d + 1)
    finite_ok = all(polyval_univar(coeffs, m) > 0 for m in range(m_min, threshold))
    tail_ok = c_d > 0 and all(c >= 0 for c in coeffs[d:])
    return {
        "name": name,
        "status": "PASS" if finite_ok and tail_ok else "FAIL",
        "dominant_degree": d,
        "dominant_coefficient": c_d,
        "lower_negative_abs_sum": neg_sum,
        "domination_threshold_integer": threshold,
        "finite_checked_range": [m_min, threshold - 1] if threshold > m_min else [],
        "finite_check_ok": finite_ok,
        "tail_coefficients_nonnegative": tail_ok,
    }


def run_telescoping_certificate() -> dict:
    b = b_poly()
    bx = deriv_x(b)
    s = s_certificate()
    r = recurrence_polynomial_r()
    lhs = add(mul(b, deriv_x(s)), mul(sub(mul(add(mpow(1), const(-3)), bx), b), s))
    rhs = mul(b, r)
    residual = sub(lhs, rhs)
    s0_ok = all(xdeg != 0 or coeff == 0 for (xdeg, _mdeg), coeff in s.items())
    return {
        "id": "all_m_telescoping_certificate",
        "status": "PASS" if not residual and s0_ok else "FAIL",
        "identity": "B*S_x + ((m-3)*B_x - B)*S = B*sum_j p_j(m)2^jB^jCbar_{m+j}",
        "B": "2*x^2 - 4*x + 1",
        "R_degrees": max_degrees(r),
        "S_degrees": max_degrees(s),
        "residual_terms": len(residual),
        "S_at_x0_zero": s0_ok,
        "boundary_infinity_zero": True,
        "valid_for": "m>=3, with m=2 checked separately",
    }


def run_edge_and_finite_sanity(max_m: int) -> dict:
    residuals = {m: recurrence_residual(m) for m in range(2, max_m + 1)}
    ok = all(value == 0 for value in residuals.values())
    return {
        "id": "edge_case_and_finite_sanity",
        "status": "PASS" if ok else "FAIL",
        "range": f"m=2..{max_m}",
        "m2_residual": residuals[2],
        "nonzero_residuals": {str(m): value for m, value in residuals.items() if value != 0},
        "first_H_values": {str(m): hm_value(m) for m in range(2, min(max_m, 10) + 1)},
    }


def run_positivity_certificate() -> dict:
    d_poly = compute_d_polynomial()
    sign_checks = [
        dominance_certificate(P0, 2, "p0_positive"),
        dominance_certificate(P1, 2, "p1_positive"),
        dominance_certificate([-c for c in P2], 2, "minus_p2_positive"),
        dominance_certificate(P3, 2, "p3_positive"),
        dominance_certificate(d_poly, 2, "D_positive"),
    ]
    base_rows = []
    base_ok = True
    for m in (2, 3, 4):
        hm = hm_value(m)
        rho = Fraction(hm_value(m + 1), hm)
        row_ok = hm > 0 and rho >= 4 * m * m
        base_ok = base_ok and row_ok
        base_rows.append(
            {
                "m": m,
                "H_m": hm,
                "rho_num": rho.numerator,
                "rho_den": rho.denominator,
                "threshold_4m2": 4 * m * m,
                "ok": row_ok,
            }
        )
    ok = base_ok and all(row["status"] == "PASS" for row in sign_checks)
    return {
        "id": "all_m_H_m_positivity_from_recurrence",
        "status": "PASS" if ok else "FAIL",
        "base_ok": base_ok,
        "base_rows": base_rows,
        "sign_and_D_checks": sign_checks,
        "D_coefficients_low_to_high": d_poly,
        "conclusion": "H_m>0 for all m>=2" if ok else "positivity certificate failed",
    }


def make_report(max_m: int) -> dict:
    checks = [
        run_telescoping_certificate(),
        run_edge_and_finite_sanity(max_m),
        run_positivity_certificate(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S8 H_m recurrence certificate",
        "status": "PASS" if not failed else "FAIL",
        "decision": "S8 H_m formula-sequence recurrence certificate passes; H_m>0 for all m>=2 is certified for the displayed formula sequence." if not failed else "S8 certificate failed; inspect failed_checks.",
        "failed_checks": failed,
        "checks": checks,
        "closed_at_S8": [
            "The order-3 recurrence for H_m is certified for all m>=3 by a concrete integration-by-parts certificate.",
            "The m=2 recurrence edge case is checked exactly.",
            "The elementary ratio induction proves H_m>0 for all m>=2 from the recurrence.",
            "S8 itself does not identify this displayed formula sequence with the true odd V_ref k=1 scalar residual; downstream S11A-H now supplies that closure.",
        ],
        "closed_downstream_after_S8": [
            "S11A-H formula-to-true-scalar closure is now closed downstream by the S11A-H edge-insertion verifier.",
        ],
        "open_after_S8": [
            "Odd V_pair-std is handled downstream by S9/S9C, and its k=1 formula-to-true-scalar link is now closed downstream by S11A-K.",
            "D~I remains Type-D scout-only; no Type-D theorem is promoted.",
        ],
        "boundary": "S8 closes the recurrence and positivity of the displayed H_m formula sequence only; downstream S11A-H owns and now supplies the formula-to-true-scalar link. S8 alone promotes no public theorem.",
    }


def print_summary(report: dict) -> None:
    print("P15-S8 H_m recurrence certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        if check["id"] == "all_m_telescoping_certificate":
            print("  residual_terms:", check["residual_terms"], "S_at_x0_zero:", check["S_at_x0_zero"])
            print("  R_degrees:", check["R_degrees"], "S_degrees:", check["S_degrees"])
        elif check["id"] == "edge_case_and_finite_sanity":
            print("  range:", check["range"], "m2_residual:", check["m2_residual"], "nonzero:", len(check["nonzero_residuals"]))
        elif check["id"] == "all_m_H_m_positivity_from_recurrence":
            print("  base_ok:", check["base_ok"])
            print("  sign/D:", [(row["name"], row["status"]) for row in check["sign_and_D_checks"]])


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S8 H_m all-m recurrence certificate")
    parser.add_argument("--max-m", type=int, default=100, help="finite sanity range after all-m certificate")
    parser.add_argument("--write-json", help="write JSON certificate")
    args = parser.parse_args()
    report = make_report(args.max_m)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
