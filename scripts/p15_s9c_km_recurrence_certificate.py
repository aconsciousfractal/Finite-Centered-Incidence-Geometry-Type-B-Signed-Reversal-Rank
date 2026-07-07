# -*- coding: utf-8 -*-
"""
P15-S9C K_m all-m recurrence and positivity certificate.

This certificate closes the final odd pair-standard k=1 residue isolated at
S9A/S9B:

    K_m = (m-1)*lambda_prime(2m+1,1) > 0,  m>=2.

For m>=4 the S9 formula is

    lambda_prime(2m+1,1)=Phi_{2m}(B(t)^(m-4)*L_m(t)),
    B(t)=2*t^2-4*t+1.

After the factorial-functional integral substitution t=1/(2x), write

    K_m = integral_0^infty exp(-x) * 2^(m+4) * A(x)^(m-4) * Q_m(x) dx,
    A(x)=2*x^2-4*x+1.

The verifier checks an explicit order-3 recurrence by an integration-by-parts
certificate over Q[m,x]. It then proves positivity by a ratio induction with
rho_m=K_{m+1}/K_m >= 4*m^2 from m>=7, plus finite base m=2..8.

No symbolic search is performed by this script. All recurrence and certificate
coefficients are fixed constants.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

SCHEMA = "p15.s9c.km_recurrence_certificate.v1"

P0 = [0, -14229504, -74274944, -148900288, -117488064, 47878144, 195155264, 202067904, 117560064, 42475520, 9469952, 1196032, 65536]
P1 = [-2327296, -5763840, 11338464, 67184080, 124526656, 125717760, 76031472, 27097664, 4785152, -8192, -143360, -16384, 0]
P2 = [-197640, -631072, 766138, 8667292, 22596688, 29952010, 20948960, 5012464, -3455168, -3383808, -1249280, -225280, -16384]
P3 = [378, 2353, -1053, -29212, -78569, -95701, -50356, 3968, 17408, 7424, 1024, 0, 0]
RECURRENCE = [P0, P1, P2, P3]
B_T_POLY = [1, -4, 2]

S_COEFFS = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [796608, 21467200, 29023568, -148022832, -453513792, -411761872, 73394032, 417847296, 302087168, 41225984, -64532992, -45281280, -13746176, -2097152, -131072],
    [-13898304, -221449728, -153372080, 2075377504, 5681716160, 5431396704, 300168864, -3548164096, -2708968496, -302898048, 659417344, 445026304, 133066752, 20086784, 1245184],
    [103386240, 1074163968, 417858016, -10957415808, -30580968096, -35307983200, -16910382656, 2080788896, 5291433344, 511321152, -1819269888, -1204641792, -355319808, -52707328, -3211264],
    [-433862784, -3460516608, -2502997856, 29453269376, 98696412736, 151613948992, 138119649728, 78712926400, 25392389600, 621631104, -4065643520, -2303664128, -664068096, -102268928, -6553600],
    [1153066752, 7764380928, 9206360896, -52059927488, -229822134784, -451010121152, -533423967168, -393467899264, -157101029248, 404263680, 38366206976, 22210396160, 6354206720, 950599680, 58982400],
    [-2064738816, -11863734528, -15036227968, 80775598400, 410574134272, 905898775872, 1175622572480, 923449089664, 376215819200, -14700691968, -108197229568, -61454467072, -17438130176, -2587492352, -159121408],
    [2529156096, 13107042304, 14768753920, -102174619648, -521944823808, -1189623119104, -1587665828608, -1269442895360, -514838360064, 35933914368, 164938236928, 92495544320, 26132545536, 3859480576, 235667456],
    [-2023252992, -10239946240, -10795458304, 84184979840, 427071573504, 983961441664, 1331806549120, 1078921136384, 439821723776, -36961341952, -149362028544, -83741868032, -23667703808, -3483762688, -210763776],
    [953607168, 4798370816, 5460402944, -37342126848, -198782620672, -474601944320, -665090198784, -557701450752, -236638666240, 16108238848, 79819358208, 45661552640, 12992249856, 1904738304, 113246208],
    [-194835456, -855616512, -1469574912, 4275642624, 35318209536, 103110117632, 166725567232, 157101958656, 75987642368, 611139584, -22484353024, -13835042816, -4028104704, -587202560, -33554432],
    [-23721984, -231723008, 56245760, 3313586176, 9026023936, 8607130112, -3066503168, -14349449728, -12649848832, -3475632128, 1759690752, 1724514304, 565182464, 82837504, 4194304],
    [18192384, 131631104, 66867712, -1438631424, -5215685632, -8664642048, -7662241280, -2881417216, 814295040, 1339662336, 544866304, 74973184, -6815744, -2097152, 0],
    [-2322432, -17166336, -11170816, 182207488, 694276096, 1210995712, 1156281344, 532568064, -32268288, -178520064, -95158272, -22544384, -2097152, 0, 0],
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


def m_shift(shift: int) -> Poly:
    return add(mpow(1), const(shift))


def max_degrees(poly: Poly) -> dict[str, int]:
    if not poly:
        return {"x": -1, "m": -1, "terms": 0}
    return {"x": max(key[0] for key in poly), "m": max(key[1] for key in poly), "terms": len(poly)}


# ---------------------------------------------------------------------------
# Univariate factorial formula for K_m.


def trim_list(poly: list[int]) -> list[int]:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_t_list_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ci in enumerate(a):
        if ci:
            for j, cj in enumerate(b):
                if cj:
                    out[i + j] += ci * cj
    return trim_list(out)


def poly_t_list_pow(poly: list[int], exponent: int) -> list[int]:
    out = [1]
    base = trim_list(poly)
    n = exponent
    while n:
        if n & 1:
            out = poly_t_list_mul(out, base)
        base = poly_t_list_mul(base, base)
        n >>= 1
    return out


def L_poly_t(m: int) -> list[int]:
    values = [0] * 9
    e = [
        2 * m - 3,
        -4 * m * m - 8 * m + 19,
        24 * m * m + 18 * m - 68,
        -56 * m * m - 24 * m + 122,
        56 * m * m + 24 * m - 92,
        -12 * m * m - 44 * m + 12,
        -16 * m * m + 56 * m,
        8 * m * m - 24 * m + 8,
    ]
    for degree, coeff in enumerate(e):
        values[degree + 1] = 2 * coeff
    return trim_list(values)


def phi_factorial(poly: list[int], n: int) -> int:
    return sum(coeff * (2 ** (n - r)) * factorial(n - r) for r, coeff in enumerate(poly) if coeff and 0 <= r <= n)


def lambda_value(m: int) -> int:
    if m == 2:
        return 28
    if m == 3:
        return 7320
    if m < 2:
        raise ValueError("lambda_value is used only for m>=2")
    poly = poly_t_list_mul(poly_t_list_pow(B_T_POLY, m - 4), L_poly_t(m))
    return phi_factorial(poly, 2 * m)


def km_value(m: int) -> int:
    return (m - 1) * lambda_value(m)


def polyval_univar(coeffs: list[int], m: int) -> int:
    total = 0
    power = 1
    for coeff in coeffs:
        total += coeff * power
        power *= m
    return total


def recurrence_residual(m: int) -> int:
    return sum(polyval_univar(poly, m) * km_value(m + j) for j, poly in enumerate(RECURRENCE))


# ---------------------------------------------------------------------------
# Integral certificate polynomials.


def e_coeff_poly(degree: int, shift: int) -> Poly:
    mm = m_shift(shift)
    mm2 = mul(mm, mm)
    rows = [
        add(scale(mm, 2), const(-3)),
        add(scale(mm2, -4), scale(mm, -8), const(19)),
        add(scale(mm2, 24), scale(mm, 18), const(-68)),
        add(scale(mm2, -56), scale(mm, -24), const(122)),
        add(scale(mm2, 56), scale(mm, 24), const(-92)),
        add(scale(mm2, -12), scale(mm, -44), const(12)),
        add(scale(mm2, -16), scale(mm, 56)),
        add(scale(mm2, 8), scale(mm, -24), const(8)),
    ]
    return rows[degree]


def qbar(shift: int) -> Poly:
    """Q_{m+shift}(x)=16*(m+shift-1)*x^8*L_{m+shift}(1/(2x))."""
    out = const(0)
    factor = add(m_shift(shift), const(-1))
    for degree in range(8):
        coeff = Fraction(2 ** (4 - degree), 1) if degree <= 4 else Fraction(1, 2 ** (degree - 4))
        term = mul(factor, e_coeff_poly(degree, shift))
        term = mul(term, xpow(7 - degree, coeff))
        out = add(out, term)
    return out


def a_poly() -> Poly:
    return add(scale(xpow(2), 2), scale(xpow(1), -4), const(1))


def s_certificate() -> Poly:
    out: Poly = {}
    for xdeg, row in enumerate(S_COEFFS):
        for mdeg, coeff in enumerate(row):
            if coeff:
                out[(xdeg, mdeg)] = Fraction(coeff, 1)
    return out


def recurrence_polynomial_r() -> Poly:
    a = a_poly()
    out = const(0)
    for shift, p_coeffs in enumerate(RECURRENCE):
        term = mul(poly_m_from_coeffs(p_coeffs), scale(pow_poly(a, shift), 2 ** shift))
        term = mul(term, qbar(shift))
        out = add(out, term)
    return out


# ---------------------------------------------------------------------------
# Positivity certificate.


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


def shift_coefficients(poly: list[int], shift: int) -> list[int]:
    """Return coefficients of poly(q+shift) in q, low-to-high."""
    out = [0] * len(poly)
    for degree, coeff in enumerate(poly):
        if coeff == 0:
            continue
        for q_degree in range(degree + 1):
            out[q_degree] += coeff * comb(degree, q_degree) * (shift ** (degree - q_degree))
    return trim_list(out)


def shifted_nonnegative_certificate(coeffs: list[int], shift: int, name: str) -> dict:
    shifted = shift_coefficients(coeffs, shift)
    all_nonnegative = all(coeff >= 0 for coeff in shifted)
    constant_positive = bool(shifted) and shifted[0] > 0
    return {
        "name": name,
        "status": "PASS" if all_nonnegative and constant_positive else "FAIL",
        "method": f"substitute m=q+{shift}; all coefficients in q are nonnegative and constant term is positive",
        "valid_for_m_ge": shift,
        "degree": len(trim_list(coeffs)) - 1,
        "shifted_degree": len(shifted) - 1,
        "shifted_constant": shifted[0] if shifted else 0,
        "shifted_coefficients_nonnegative": all_nonnegative,
        "shifted_constant_positive": constant_positive,
        "negative_shifted_coefficients": {str(i): coeff for i, coeff in enumerate(shifted) if coeff < 0},
    }


def compute_ratio_d_polynomial() -> list[int]:
    t_m = [0, 0, 4]
    t_m1 = [4, 8, 4]
    t_m2 = [16, 16, 4]
    term_main = mul_list(mul_list(t_m, t_m1), scale_list(P2, -1))
    term_p0 = P0
    term_p3 = mul_list(mul_list(mul_list(t_m, t_m1), t_m2), P3)
    return add_list(add_list(term_main, term_p0, -1), term_p3, -1)


def run_telescoping_certificate() -> dict:
    a = a_poly()
    ax = deriv_x(a)
    s = s_certificate()
    r = recurrence_polynomial_r()
    lhs = add(mul(a, deriv_x(s)), mul(sub(mul(add(mpow(1), const(-4)), ax), a), s))
    rhs = mul(a, r)
    residual = sub(lhs, rhs)
    s0_ok = all(xdeg != 0 or coeff == 0 for (xdeg, _mdeg), coeff in s.items())
    return {
        "id": "all_m_telescoping_certificate",
        "status": "PASS" if not residual and s0_ok else "FAIL",
        "identity": "A*S_x + ((m-4)*A_x - A)*S = A*sum_j p_j(m)2^jA^jQ_{m+j}",
        "A": "2*x^2 - 4*x + 1",
        "R_degrees": max_degrees(r),
        "S_degrees": max_degrees(s),
        "residual_terms": len(residual),
        "S_at_x0_zero": s0_ok,
        "boundary_infinity_zero": True,
        "valid_for": "m>=4, with m=2..3 recurrence edge cases checked exactly",
    }


def run_edge_and_finite_sanity(max_m: int) -> dict:
    residuals = {m: recurrence_residual(m) for m in range(2, max_m + 1)}
    ok = all(value == 0 for value in residuals.values())
    return {
        "id": "edge_case_and_finite_sanity",
        "status": "PASS" if ok else "FAIL",
        "range": f"m=2..{max_m}",
        "nonzero_residuals": {str(m): value for m, value in residuals.items() if value != 0},
        "first_K_values": {str(m): km_value(m) for m in range(2, min(max_m, 10) + 1)},
        "first_lambda_values": {str(m): lambda_value(m) for m in range(2, min(max_m, 10) + 1)},
    }


def run_positivity_certificate() -> dict:
    d_poly = compute_ratio_d_polynomial()
    sign_checks = [
        shifted_nonnegative_certificate(P0, 7, "p0_positive_from_7"),
        shifted_nonnegative_certificate(scale_list(P1, -1), 7, "p1_nonpositive_from_7"),
        shifted_nonnegative_certificate(scale_list(P2, -1), 7, "minus_p2_positive_from_7"),
        shifted_nonnegative_certificate(P3, 7, "p3_positive_from_7"),
        shifted_nonnegative_certificate(d_poly, 7, "ratio_D_positive_from_7"),
    ]
    base_rows = []
    base_ok = True
    for m in range(2, 9):
        km = km_value(m)
        rho = Fraction(km_value(m + 1), km)
        threshold = 4 * m * m
        row_ok = km > 0 and rho >= threshold
        base_ok = base_ok and row_ok
        base_rows.append(
            {
                "m": m,
                "K_m": km,
                "rho_num": rho.numerator,
                "rho_den": rho.denominator,
                "threshold_4m2": threshold,
                "ok": row_ok,
            }
        )
    ok = base_ok and all(row["status"] == "PASS" for row in sign_checks)
    return {
        "id": "all_m_K_m_positivity_from_recurrence",
        "status": "PASS" if ok else "FAIL",
        "base_range": "m=2..8 for K_m>0 and rho_m>=4m^2",
        "base_ok": base_ok,
        "base_rows": base_rows,
        "sign_and_D_checks": sign_checks,
        "D_coefficients_low_to_high": d_poly,
        "ratio_induction_start_m": 7,
        "ratio_threshold": "rho_m>=4*m^2",
        "conclusion": "K_m>0 for all m>=2" if ok else "positivity certificate failed",
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
        "gate": "P15-S9C K_m recurrence certificate",
        "status": "PASS" if not failed else "FAIL",
        "decision": "S9C K_m formula-sequence recurrence certificate passes; K_m>0 for all m>=2 is certified for the displayed formula sequence." if not failed else "S9C certificate failed; inspect failed_checks.",
        "failed_checks": failed,
        "checks": checks,
        "closed_at_S9C": [
            "The order-3 recurrence for K_m is certified for all m>=4 by a concrete integration-by-parts certificate.",
            "The recurrence edge cases m=2..3 are checked exactly by finite sanity.",
            "The ratio induction proves K_m>0 for all m>=2 from the recurrence and finite base rows.",
            "S9C itself does not identify this displayed formula sequence with the true odd pair-standard k=1 scalar; downstream S11A-K now supplies that closure.",
        ],
        "closed_downstream_after_S9C": [
            "S11A-K formula-to-true-scalar closure is now closed downstream by the S11A-K global product/Phi certificate.",
        ],
        "open_after_S9C": [
            "Positive-only separator and full B_n fingerprint remain separate research tasks.",
            "D~I remains Type-D scout-only; no Type-D theorem is promoted.",
        ],
        "boundary": "S9C closes recurrence and positivity of the displayed K_m formula sequence only; downstream S11A-K owns and now supplies the formula-to-true-scalar link. S9C alone promotes no public theorem or Type-D/fingerprint claim.",
    }


def print_summary(report: dict) -> None:
    print("P15-S9C K_m recurrence certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        if check["id"] == "all_m_telescoping_certificate":
            print("  residual_terms:", check["residual_terms"], "S_at_x0_zero:", check["S_at_x0_zero"])
            print("  R_degrees:", check["R_degrees"], "S_degrees:", check["S_degrees"])
        elif check["id"] == "edge_case_and_finite_sanity":
            print("  range:", check["range"], "nonzero:", len(check["nonzero_residuals"]))
            first = check["first_K_values"]
            print("  K_2..K_5:", [first[str(m)] for m in range(2, 6)])
        elif check["id"] == "all_m_K_m_positivity_from_recurrence":
            print("  base_ok:", check["base_ok"], "start_m:", check["ratio_induction_start_m"])
            print("  sign/D:", [(row["name"], row["status"]) for row in check["sign_and_D_checks"]])


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S9C K_m all-m recurrence certificate")
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
