"""P15-S11A symbolic transfer verifier.

This is intentionally not a search and not a held-out replay. It verifies the
local polynomial algebra that turns the S11A intermediate coefficient
extractions into the displayed univariate factorial formula targets.

Boundary: this script is a local transfer certificate. The all-m true-scalar
links are closed downstream by the S11A-H edge-insertion certificate and the
S11A-K global product/Phi certificate.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

SCHEMA = "p15.s11a.scalar_link_symbolic_verifier.v1"

# Polynomial ring Z[m,t], represented by {(m_degree, t_degree): coefficient}.
Key = tuple[int, int]
Poly = dict[Key, int]
Comp = dict[str, Poly]


def clean(poly: Poly) -> Poly:
    return {key: value for key, value in poly.items() if value}


def const(value: int) -> Poly:
    return {} if value == 0 else {(0, 0): value}


def m_poly(coeff: int = 1, shift: int = 0) -> Poly:
    return clean({(1, 0): coeff, (0, 0): shift})


def t_poly(coeffs: list[int]) -> Poly:
    return clean({(0, degree): coeff for degree, coeff in enumerate(coeffs) if coeff})


def add(*polys: Poly) -> Poly:
    out: dict[Key, int] = defaultdict(int)
    for poly in polys:
        for key, value in poly.items():
            out[key] += value
    return clean(dict(out))


def neg(poly: Poly) -> Poly:
    return {key: -value for key, value in poly.items()}


def sub(a: Poly, b: Poly) -> Poly:
    return add(a, neg(b))


def scale(poly: Poly, scalar: int) -> Poly:
    return clean({key: scalar * value for key, value in poly.items()})


def mul(a: Poly, b: Poly) -> Poly:
    out: dict[Key, int] = defaultdict(int)
    for (am, at), av in a.items():
        for (bm, bt), bv in b.items():
            out[(am + bm, at + bt)] += av * bv
    return clean(dict(out))


def pow_small(poly: Poly, exponent: int) -> Poly:
    out = const(1)
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def terms(poly: Poly) -> list[dict[str, int]]:
    return [
        {"m_degree": md, "t_degree": td, "coefficient": coeff}
        for (md, td), coeff in sorted(poly.items(), key=lambda item: (item[0][1], item[0][0]))
    ]


def max_degrees(poly: Poly) -> dict[str, int]:
    if not poly:
        return {"m": -1, "t": -1, "terms": 0}
    return {"m": max(md for md, _ in poly), "t": max(td for _, td in poly), "terms": len(poly)}


ZERO = const(0)
ONE = const(1)
M = m_poly()
T = t_poly([0, 1])
A = t_poly([1, -4, 2])
B = A
C = t_poly([0, 2, -2])  # 2t-2t^2: coefficient of u and v in the free-pair transfer.
ONE_MINUS_T = t_poly([1, -1])
ONE_MINUS_2T = t_poly([1, -2])
TWO_T_MINUS_ONE = t_poly([-1, 2])


def comp(constant: Poly | None = None, u: Poly | None = None, v: Poly | None = None, uv: Poly | None = None) -> Comp:
    return {
        "const": constant or ZERO,
        "u": u or ZERO,
        "v": v or ZERO,
        "uv": uv or ZERO,
    }


def comp_add(*items: Comp) -> Comp:
    return {name: add(*(item[name] for item in items)) for name in ("const", "u", "v", "uv")}


def comp_scale(item: Comp, scalar: int) -> Comp:
    return {name: scale(item[name], scalar) for name in item}


def comp_sub(a: Comp, b: Comp) -> Comp:
    return comp_add(a, comp_scale(b, -1))


def comp_mul(a: Comp, b: Comp) -> Comp:
    # Multiplication in Z[m,t][u,v]/(u^2,v^2), enough for [uv].
    return comp(
        constant=mul(a["const"], b["const"]),
        u=add(mul(a["u"], b["const"]), mul(a["const"], b["u"])),
        v=add(mul(a["v"], b["const"]), mul(a["const"], b["v"])),
        uv=add(
            mul(a["uv"], b["const"]),
            mul(a["u"], b["v"]),
            mul(a["v"], b["u"]),
            mul(a["const"], b["uv"]),
        ),
    )


def h_cubic() -> Poly:
    # 2m t^3 -(4m+2)t^2 +(2m+2)t - 1.
    return add(
        const(-1),
        mul(add(scale(M, 2), const(2)), T),
        mul(add(scale(M, -4), const(-2)), pow_small(T, 2)),
        mul(scale(M, 2), pow_small(T, 3)),
    )


def h_transfer_identity() -> dict:
    # Agent-1 intermediate:
    # R_m(z,t)=az(z,t)^(m-2) (1+(z-2)t)(1-t)(1-2t),
    # S_m=Phi(R_m(1,t)-d_z R_m(1,t)).  Factoring A(t)^(m-3),
    # this local polynomial should equal (1-t)(2t-1)*cubic.
    m_minus_2 = sub(M, const(2))
    dz_a_at_1 = C
    tail_at_1 = ONE_MINUS_T
    dz_tail = T
    bracket = sub(sub(mul(A, tail_at_1), mul(mul(m_minus_2, dz_a_at_1), tail_at_1)), mul(A, dz_tail))
    local_from_agent1 = mul(mul(ONE_MINUS_T, ONE_MINUS_2T), bracket)
    target_l = mul(mul(ONE_MINUS_T, TWO_T_MINUS_ONE), h_cubic())
    residual_l = sub(local_from_agent1, target_l)

    p_m = sub(mul(scale(sub(M, const(1)), 2), mul(TWO_T_MINUS_ONE, h_cubic())), pow_small(A, 2))
    local_h = sub(mul(scale(sub(M, const(1)), 2), target_l), mul(ONE_MINUS_T, pow_small(A, 2)))
    target_h = mul(ONE_MINUS_T, p_m)
    residual_h = sub(local_h, target_h)
    return {
        "id": "H_transfer_agent1_to_univariate_P_m",
        "status": "PASS" if not residual_l and not residual_h else "FAIL",
        "description": "Verifies R_m(1,t)-d_z R_m(1,t) local algebra and the resulting P_m target over Z[m,t].",
        "residual_L_terms": terms(residual_l),
        "residual_H_terms": terms(residual_h),
        "target_L_degrees": max_degrees(target_l),
        "target_H_degrees": max_degrees(target_h),
    }


def k_l_poly() -> Poly:
    rows = [
        add(scale(M, 2), const(-3)),
        add(scale(mul(M, M), -4), scale(M, -8), const(19)),
        add(scale(mul(M, M), 24), scale(M, 18), const(-68)),
        add(scale(mul(M, M), -56), scale(M, -24), const(122)),
        add(scale(mul(M, M), 56), scale(M, 24), const(-92)),
        add(scale(mul(M, M), -12), scale(M, -44), const(12)),
        add(scale(mul(M, M), -16), scale(M, 56)),
        add(scale(mul(M, M), 8), scale(M, -24), const(8)),
    ]
    out = ZERO
    for i, row in enumerate(rows, start=1):
        out = add(out, mul(scale(row, 2), pow_small(T, i)))
    return out


def free_pair_power_components(exponent_shift: int) -> Comp:
    # For base b+c(u+v), compute components of base^(m+shift), after
    # factoring b^(m-4).  Only shifts -1 and -2 are needed here.
    if exponent_shift == -1:
        n_poly = sub(M, const(1))
        n_minus_1 = sub(M, const(2))
        return comp(
            constant=pow_small(B, 3),
            u=mul(mul(n_poly, pow_small(B, 2)), C),
            v=mul(mul(n_poly, pow_small(B, 2)), C),
            uv=mul(mul(mul(n_poly, n_minus_1), B), pow_small(C, 2)),
        )
    if exponent_shift == -2:
        n_poly = sub(M, const(2))
        n_minus_1 = sub(M, const(3))
        return comp(
            constant=pow_small(B, 2),
            u=mul(mul(n_poly, B), C),
            v=mul(mul(n_poly, B), C),
            uv=mul(mul(n_poly, n_minus_1), pow_small(C, 2)),
        )
    raise ValueError("unsupported exponent shift")


def k_transfer_identity() -> dict:
    # Starting from the S9 bivariate target:
    # c * ( a^(m-1)*((1+u)sx+(1+v)sy) - 4*a^(m-2)*h^2 ).
    # After extracting [uv] and factoring B(t)^(m-4), the local polynomial
    # must be exactly the displayed degree-eight L_m(t).
    a_m1 = free_pair_power_components(-1)
    a_m2 = free_pair_power_components(-2)

    d_endpoint = comp(constant=scale(ONE_MINUS_T, 2), u=ONE, v=ONE)
    h = comp(constant=t_poly([1, -2]), u=T, v=T)
    h_squared = comp_mul(h, h)
    bracket = comp_sub(comp_mul(a_m1, d_endpoint), comp_scale(comp_mul(a_m2, h_squared), 4))
    c_endpoint = comp(constant=ONE_MINUS_T, uv=T)
    expr = comp_mul(c_endpoint, bracket)
    local = expr["uv"]
    target = k_l_poly()
    residual = sub(local, target)
    return {
        "id": "K_bivariate_to_univariate_L_m",
        "status": "PASS" if not residual else "FAIL",
        "description": "Verifies the S9 bivariate k=1 target collapses to B(t)^(m-4)L_m(t) over Z[m,t].",
        "residual_terms": terms(residual),
        "local_degrees": max_degrees(local),
        "target_degrees": max_degrees(target),
    }


def make_report() -> dict:
    checks = [h_transfer_identity(), k_transfer_identity()]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S11A scalar-link symbolic transfer verifier",
        "status": "PASS" if not failed else "FAIL",
        "decision": (
            "Symbolic transfer identities pass over Z[m,t]; S11A formula targets are algebraically coherent. "
            "Downstream S11A-H and S11A-K certificates now close both true-scalar all-m links."
            if not failed
            else "Symbolic transfer verifier failed; inspect residuals."
        ),
        "failed_checks": failed,
        "checks": checks,
        "no_go_rules_obeyed": {
            "no_held_out_extension": True,
            "no_symbolic_search": True,
            "does_not_use_formula_sequence_as_true_scalar_check": True,
        },
        "closed_by_this_script": [
            "H_m intermediate transfer algebra to the displayed P_m polynomial.",
            "K_m S9 bivariate target transfer algebra to the displayed L_m polynomial.",
        ],
        "closed_downstream": [
            "S11A-H is closed by P15_S11A_H_TRUE_SCALAR_EDGE_INSERTION_VERIFIER plus S8.",
            "S11A-K is closed by P15_S11A_K_TRUE_SCALAR_BIVARIATE_VERIFIER plus P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE plus S9C.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S11A scalar-link symbolic transfer verifier")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        for key in ("target_L_degrees", "target_H_degrees", "local_degrees", "target_degrees"):
            if key in check:
                print(" ", key + ":", check[key])
        residual_count = sum(len(check.get(name, [])) for name in ("residual_terms", "residual_L_terms", "residual_H_terms"))
        print("  residual_terms:", residual_count)


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S11A symbolic transfer verifier")
    parser.add_argument("--write-json", help="write JSON report")
    args = parser.parse_args()
    report = make_report()
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
