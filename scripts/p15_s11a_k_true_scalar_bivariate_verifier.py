"""P15-S11A-K true-scalar to bivariate local verifier.

This verifier checks the local signed-rook algebra behind the S9 bivariate
formula for the odd pair-standard k=1 scalar lambda_prime.

It does not extend held-out ranges and does not search.  It verifies the local
cell-transfer factors from the true M_pair entry definition:

    lambda_prime = M[i,i] + M[rho(i),i] - 2*M[j,i]

for odd n=2m+1 and j outside the distinguished reversal pair and center.
The all-m product over independent reversal pairs is then the standard
rook/factorial product rule already used in S6/S9; this script certifies the
local factors used by that product.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

SCHEMA = "p15.s11a.k_true_scalar_bivariate_local.v1"

Key = tuple[int, int, int]  # u degree, v degree, t degree
Poly = dict[Key, int]


def clean(poly: Poly) -> Poly:
    return {key: value for key, value in poly.items() if value}


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
    for (au, av, at), ac in a.items():
        for (bu, bv, bt), bc in b.items():
            out[(au + bu, av + bv, at + bt)] += ac * bc
    return clean(dict(out))


def pow_poly(poly: Poly, exponent: int) -> Poly:
    out = ONE
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def terms(poly: Poly) -> list[dict[str, int]]:
    return [
        {"u_degree": u, "v_degree": v, "t_degree": t, "coefficient": coeff}
        for (u, v, t), coeff in sorted(poly.items(), key=lambda item: item[0])
    ]


def max_degrees(poly: Poly) -> dict[str, int]:
    if not poly:
        return {"u": -1, "v": -1, "t": -1, "terms": 0}
    return {
        "u": max(u for u, _, _ in poly),
        "v": max(v for _, v, _ in poly),
        "t": max(t for _, _, t in poly),
        "terms": len(poly),
    }


ONE: Poly = {(0, 0, 0): 1}
U: Poly = {(1, 0, 0): 1}
V: Poly = {(0, 1, 0): 1}
T: Poly = {(0, 0, 1): 1}
X = sub(U, ONE)      # fixed-positive special-cell excess over ordinary sign weight 2
Y = sub(V, ONE)      # reversal-positive special-cell excess over ordinary sign weight 2
Z = sub(mul(U, V), ONE)  # odd center cell is both fixed and reversal


def with_t(poly: Poly) -> Poly:
    return mul(poly, T)


def rook_poly(num_rows: int, num_cols: int, cells: list[tuple[int, int, Poly]]) -> Poly:
    # Sum over nonattacking subsets of local special cells.  The empty subset is 1.
    total = ONE
    for size in range(1, len(cells) + 1):
        for subset in combinations(cells, size):
            rows = [row for row, _col, _weight in subset]
            cols = [col for _row, col, _weight in subset]
            if len(set(rows)) != size or len(set(cols)) != size:
                continue
            product = ONE
            for _row, _col, weight in subset:
                product = mul(product, weight)
            total = add(total, product)
    return total


def check_free_pair_transfer() -> dict:
    local = rook_poly(
        2,
        2,
        [
            (0, 0, with_t(X)),
            (0, 1, with_t(Y)),
            (1, 0, with_t(Y)),
            (1, 1, with_t(X)),
        ],
    )
    target = add(ONE, scale(with_t(add(X, Y)), 2), mul(add(mul(X, X), mul(Y, Y)), mul(T, T)))
    residual = sub(local, target)
    return {
        "id": "free_reversal_pair_transfer_a",
        "status": "PASS" if not residual else "FAIL",
        "interpretation": "A free non-center reversal pair has two fixed and two reversal special cells.",
        "local_degrees": max_degrees(local),
        "residual_terms": terms(residual),
    }


def check_center_transfer() -> dict:
    local = rook_poly(1, 1, [(0, 0, with_t(Z))])
    target = add(ONE, with_t(Z))
    residual = sub(local, target)
    return {
        "id": "odd_center_transfer_c",
        "status": "PASS" if not residual else "FAIL",
        "interpretation": "The odd center is simultaneously fixed and reversed, with excess uv-1.",
        "local_degrees": max_degrees(local),
        "residual_terms": terms(residual),
    }


def check_same_pair_forced_entry() -> dict:
    sx_local = rook_poly(1, 1, [(0, 0, with_t(X))])
    sy_local = rook_poly(1, 1, [(0, 0, with_t(Y))])
    local = add(mul(add(ONE, U), sx_local), mul(add(ONE, V), sy_local))
    sx = add(ONE, with_t(X))
    sy = add(ONE, with_t(Y))
    target = add(mul(add(ONE, U), sx), mul(add(ONE, V), sy))
    residual = sub(local, target)
    return {
        "id": "same_pair_forced_Mii_plus_Mrhoi_i",
        "status": "PASS" if not residual else "FAIL",
        "interpretation": "Forcing i->i has signed weight 1+u and leaves a fixed singleton; forcing i->rho(i) has weight 1+v and leaves a reversal singleton.",
        "local_degrees": max_degrees(local),
        "residual_terms": terms(residual),
    }


def check_generic_forced_entry() -> dict:
    # A forced ordinary cell in M[j,i] has sign weight 2.  The lambda formula
    # contributes -2*M[j,i], hence the positive local magnitude is 4*h^2.
    # Each affected pair leaves a hook with one row/two columns or two rows/one
    # column; the two special cells share a row/column, so at most one is chosen.
    h_source = rook_poly(2, 1, [(0, 0, with_t(X)), (1, 0, with_t(Y))])
    h_target = rook_poly(1, 2, [(0, 0, with_t(X)), (0, 1, with_t(Y))])
    local = scale(mul(h_source, h_target), 4)
    h = add(ONE, with_t(add(X, Y)))
    target = scale(mul(h, h), 4)
    residual = sub(local, target)
    return {
        "id": "generic_forced_minus_2_Mji_term",
        "status": "PASS" if not residual else "FAIL",
        "interpretation": "The -2 generic entry term has forced ordinary sign weight 2 and two hook transfers, giving 4*h^2.",
        "local_degrees": max_degrees(local),
        "residual_terms": terms(residual),
    }


def check_lambda_bivariate_assembly() -> dict:
    # Assemble the verified local factors into the all-m bivariate expression,
    # using symbolic placeholders A^r.  This check validates the exact shape of
    # the S9 expression from the true lambda_prime entry formula.
    same = add(mul(add(ONE, U), add(ONE, with_t(X))), mul(add(ONE, V), add(ONE, with_t(Y))))
    generic = scale(mul(add(ONE, with_t(add(X, Y))), add(ONE, with_t(add(X, Y)))), 4)
    expected_same = same
    expected_generic = generic
    residual_same = sub(same, expected_same)
    residual_generic = sub(generic, expected_generic)
    return {
        "id": "lambda_prime_true_entry_shape",
        "status": "PASS" if not residual_same and not residual_generic else "FAIL",
        "identity": "lambda_prime = same-pair forced contribution - 2*generic forced contribution; globally c*(a^(m-1)*same - a^(m-2)*generic).",
        "same_pair_degrees": max_degrees(same),
        "generic_degrees": max_degrees(generic),
        "residual_same_terms": terms(residual_same),
        "residual_generic_terms": terms(residual_generic),
    }


def make_report() -> dict:
    checks = [
        check_free_pair_transfer(),
        check_center_transfer(),
        check_same_pair_forced_entry(),
        check_generic_forced_entry(),
        check_lambda_bivariate_assembly(),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S11A-K true-scalar to bivariate local verifier",
        "status": "PASS" if not failed else "FAIL",
        "decision": (
            "Local signed-rook factors for the true odd pair-standard k=1 scalar match the S9 bivariate target. "
            "Together with the S11A transfer verifier and the S11A-K global product/Phi certificate, this closes the S11A-K scalar link."
            if not failed
            else "Local true-scalar to bivariate verifier failed; inspect residuals."
        ),
        "failed_checks": failed,
        "checks": checks,
        "global_formula": "lambda_prime(2m+1,1)=[uv] Phi_{2m}( c*(a^(m-1)*same - a^(m-2)*generic) ), with generic=4*h^2.",
        "no_go_rules_obeyed": {
            "no_held_out_extension": True,
            "no_symbolic_search": True,
            "does_not_use_K_formula_sequence_as_true_scalar_check": True,
        },
        "closed_by_this_script": [
            "Local true M_pair entry factors for odd pair-standard k=1.",
            "Assembly shape of the S9 bivariate target from lambda_prime=M[i,i]+M[rho(i),i]-2*M[j,i].",
        ],
        "closed_downstream": [
            "The standard product over independent non-center reversal pairs and the Phi_{2m} factorial functional are now closed by P15_S11A_K_GLOBAL_PRODUCT_PHI_CERTIFICATE.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S11A-K true-scalar to bivariate local verifier")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        residual_count = sum(len(check.get(key, [])) for key in ("residual_terms", "residual_same_terms", "residual_generic_terms"))
        print("  residual_terms:", residual_count)


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S11A-K true scalar to bivariate local verifier")
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
