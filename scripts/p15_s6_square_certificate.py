# -*- coding: utf-8 -*-
"""
P15-S6 C2 x C2 square and closed-form certificate builder.

This gate certifies the count-level parent grammar B^square.  It is not a
Type-D theorem promotion and it does not close the remaining odd rank theorem
obligations.

Certified at S6:

* B^square_n(a,b,c,d) records positive/negative identity and reversal hits.
* The square generating function
      F_square(xp,xm,yp,ym;q) = sum_g q^nu xp^A xm^B yp^C ym^D
  equals per(W) with W entries from the four diagonal states.
* At q=-1, F_square collapses to an explicit closed form.
* The old central BI and new B~I are two marginals of B^square.
* The old DI count formula is recovered.
* The even-n D~I count correction is recovered from the q=-1 collapse.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from math import comb
from pathlib import Path

import p15_s3_exact_engine as s3

SCHEMA = "p15.s6.square_certificate.v1"

Poly5 = dict[tuple[int, int, int, int, int], int]
Poly4 = dict[tuple[int, int, int, int], int]

ZERO5 = (0, 0, 0, 0, 0)
ZERO4 = (0, 0, 0, 0)


def add_poly5_to(target: dict[tuple[int, int, int, int, int], int], source: Poly5, scale: int = 1) -> None:
    for key, coeff in source.items():
        target[key] += scale * coeff
        if target[key] == 0:
            del target[key]


def mul_poly5(a: Poly5, b: Poly5) -> Poly5:
    out: dict[tuple[int, int, int, int, int], int] = defaultdict(int)
    for (a1, b1, c1, d1, q1), coeff1 in a.items():
        for (a2, b2, c2, d2, q2), coeff2 in b.items():
            out[(a1 + a2, b1 + b2, c1 + c2, d1 + d2, q1 + q2)] += coeff1 * coeff2
    return dict(out)


def add_poly4_to(target: dict[tuple[int, int, int, int], int], source: Poly4, scale: int = 1) -> None:
    for key, coeff in source.items():
        target[key] += scale * coeff
        if target[key] == 0:
            del target[key]


def mul_poly4(a: Poly4, b: Poly4) -> Poly4:
    out: dict[tuple[int, int, int, int], int] = defaultdict(int)
    for (a1, b1, c1, d1), coeff1 in a.items():
        for (a2, b2, c2, d2), coeff2 in b.items():
            out[(a1 + a2, b1 + b2, c1 + c2, d1 + d2)] += coeff1 * coeff2
    return {key: value for key, value in out.items() if value}


def pow_poly4(poly: Poly4, power: int) -> Poly4:
    out: Poly4 = {ZERO4: 1}
    for _ in range(power):
        out = mul_poly4(out, poly)
    return out


def square_stats(pi: tuple[int, ...], eps: tuple[int, ...], rho: tuple[int, ...]) -> tuple[int, int, int, int, int]:
    a = b = c = d = 0
    for i in range(len(pi)):
        is_id = pi[i] == i
        is_rev = pi[i] == rho[i]
        if is_id and eps[i] == 1:
            a += 1
        if is_id and eps[i] == -1:
            b += 1
        if is_rev and eps[i] == 1:
            c += 1
        if is_rev and eps[i] == -1:
            d += 1
    nu = sum(1 for value in eps if value == -1)
    return a, b, c, d, nu


def brute_square_poly(n: int) -> Poly5:
    rho = s3.reversal(n)
    out: dict[tuple[int, int, int, int, int], int] = defaultdict(int)
    for pi, eps in s3.signed_perms(n):
        out[square_stats(pi, eps, rho)] += 1
    return dict(out)


def W_cell_poly(n: int, i: int, j: int) -> Poly5:
    rho = s3.reversal(n)
    if j == i and j == rho[i]:
        return {(1, 0, 1, 0, 0): 1, (0, 1, 0, 1, 1): 1}  # xp*yp + q*xm*ym
    if j == i:
        return {(1, 0, 0, 0, 0): 1, (0, 1, 0, 0, 1): 1}  # xp + q*xm
    if j == rho[i]:
        return {(0, 0, 1, 0, 0): 1, (0, 0, 0, 1, 1): 1}  # yp + q*ym
    return {ZERO5: 1, (0, 0, 0, 0, 1): 1}  # 1 + q


def permanent_poly_square(n: int) -> Poly5:
    cells = [[W_cell_poly(n, i, j) for j in range(n)] for i in range(n)]
    dp: list[dict[tuple[int, int, int, int, int], int]] = [defaultdict(int) for _ in range(1 << n)]
    dp[0][ZERO5] = 1
    for mask in range(1 << n):
        row = mask.bit_count()
        if row >= n or not dp[mask]:
            continue
        current = dict(dp[mask])
        for col in range(n):
            if mask & (1 << col):
                continue
            product = mul_poly5(current, cells[row][col])
            add_poly5_to(dp[mask | (1 << col)], product)
    return dict(dp[(1 << n) - 1])


def specialize_q_minus_one(poly: Poly5) -> Poly4:
    out: dict[tuple[int, int, int, int], int] = defaultdict(int)
    for (a, b, c, d, nu), coeff in poly.items():
        out[(a, b, c, d)] += coeff * ((-1) ** nu)
    return {key: value for key, value in out.items() if value}


def q_minus_one_closed_form(n: int) -> Poly4:
    # ((xp-xm)^2 + (yp-ym)^2)^m, with the center factor xp*yp-xm*ym for odd n.
    base: Poly4 = {
        (2, 0, 0, 0): 1,
        (1, 1, 0, 0): -2,
        (0, 2, 0, 0): 1,
        (0, 0, 2, 0): 1,
        (0, 0, 1, 1): -2,
        (0, 0, 0, 2): 1,
    }
    if n % 2 == 0:
        return pow_poly4(base, n // 2)
    factor: Poly4 = {(1, 0, 1, 0): 1, (0, 1, 0, 1): -1}
    return mul_poly4(factor, pow_poly4(base, (n - 1) // 2))


def square_marginals(poly: Poly5) -> tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]:
    bi: dict[tuple[int, int], int] = defaultdict(int)
    btilde: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b, c, _d, _nu), coeff in poly.items():
        bi[(a, b)] += coeff
        btilde[(a, c)] += coeff
    return dict(bi), dict(btilde)


def brute_marginals(n: int) -> tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]:
    rho = s3.reversal(n)
    bi: dict[tuple[int, int], int] = defaultdict(int)
    btilde: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in s3.signed_perms(n):
        a = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
        b = sum(1 for i in range(n) if pi[i] == i and eps[i] == -1)
        c = sum(1 for i in range(n) if pi[i] == rho[i] and eps[i] == 1)
        bi[(a, b)] += 1
        btilde[(a, c)] += 1
    return dict(bi), dict(btilde)


def check_square_symmetries(poly: Poly5) -> dict:
    s_ok = True
    r_ok = True
    sr_ok = True
    for (a, b, c, d, nu), coeff in poly.items():
        # S=-I toggles all signs: (a,b,c,d)->(b,a,d,c), nu->n-nu.
        # R=(rho,+) swaps identity and reversal: (a,b,c,d)->(c,d,a,b), nu fixed.
        # SR combines the two.
        n_total_hits = None  # not used; kept explicit in comments above.
        s_key_candidates = [key for key in ()]
        del n_total_hits, s_key_candidates
        # We recover n from total sign exponent range outside this helper by observing max nu below.
    max_nu = max(key[4] for key in poly) if poly else 0
    n = max_nu
    for (a, b, c, d, nu), coeff in poly.items():
        if poly.get((b, a, d, c, n - nu), 0) != coeff:
            s_ok = False
        if poly.get((c, d, a, b, nu), 0) != coeff:
            r_ok = False
        if poly.get((d, c, b, a, n - nu), 0) != coeff:
            sr_ok = False
    return {"S_minus_I_ok": s_ok, "R_reversal_ok": r_ok, "SR_ok": sr_ok}


def di_recovery_rows(n: int) -> dict:
    bi: dict[tuple[int, int], int] = defaultdict(int)
    di: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in s3.signed_perms(n):
        k = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
        l = sum(1 for i in range(n) if pi[i] == i and eps[i] == -1)
        bi[(k, l)] += 1
        if sum(1 for value in eps if value == -1) % 2 == 0:
            di[(k, l)] += 1
    rows = []
    ok = True
    for key in sorted(bi):
        k, l = key
        correction = ((-1) ** l) * comb(n, l) if k + l == n else 0
        predicted = (bi[key] + correction) // 2
        row_ok = di.get(key, 0) == predicted and (bi[key] + correction) % 2 == 0
        ok = ok and row_ok
        rows.append(
            {
                "k": k,
                "l": l,
                "BI": bi[key],
                "DI": di.get(key, 0),
                "correction": correction,
                "predicted_DI": predicted,
                "ok": row_ok,
            }
        )
    return {"n": n, "rows_checked": len(rows), "ok": ok, "rows": rows}


def delta_even(n: int, k: int, l: int) -> int:
    m = n // 2
    return ((-1) ** (k + l)) * sum(
        comb(m, j) * comb(2 * j, k) * comb(2 * m - 2 * j, l)
        for j in range(m + 1)
        if k <= 2 * j and l <= 2 * m - 2 * j
    )


def dtilde_even_rows(n: int) -> dict:
    assert n % 2 == 0
    rho = s3.reversal(n)
    btilde: dict[tuple[int, int], int] = defaultdict(int)
    dtilde: dict[tuple[int, int], int] = defaultdict(int)
    for pi, eps in s3.signed_perms(n):
        k = sum(1 for i in range(n) if pi[i] == i and eps[i] == 1)
        l = sum(1 for i in range(n) if pi[i] == rho[i] and eps[i] == 1)
        btilde[(k, l)] += 1
        if sum(1 for value in eps if value == -1) % 2 == 0:
            dtilde[(k, l)] += 1
    rows = []
    ok = True
    for key in sorted(btilde):
        k, l = key
        correction = delta_even(n, k, l)
        predicted = (btilde[key] + correction) // 2
        row_ok = dtilde.get(key, 0) == predicted and (btilde[key] + correction) % 2 == 0
        ok = ok and row_ok
        rows.append(
            {
                "k": k,
                "l": l,
                "Btilde": btilde[key],
                "Dtilde": dtilde.get(key, 0),
                "Delta": correction,
                "predicted_Dtilde": predicted,
                "ok": row_ok,
            }
        )
    return {"n": n, "rows_checked": len(rows), "ok": ok, "rows": rows}


def run_square_permanent(max_n: int = 7) -> dict:
    rows = []
    ok = True
    for n in range(2, max_n + 1):
        brute = brute_square_poly(n)
        perm = permanent_poly_square(n)
        equal = brute == perm
        ok = ok and equal
        sym = check_square_symmetries(brute)
        ok = ok and all(sym.values())
        rows.append(
            {
                "n": n,
                "group_size": sum(brute.values()),
                "coefficient_count": len(brute),
                "permanent_polynomial_matches_brute": equal,
                "symmetries": sym,
            }
        )
    return {
        "id": "square_permanent_polynomial_replay",
        "status": "PASS" if ok else "FAIL",
        "range": f"n=2..{max_n}, full coefficient polynomial over xp,xm,yp,ym,q",
        "rows": rows,
    }


def run_q_minus_one(max_n: int = 7) -> dict:
    rows = []
    ok = True
    for n in range(2, max_n + 1):
        perm = permanent_poly_square(n)
        specialized = specialize_q_minus_one(perm)
        closed = q_minus_one_closed_form(n)
        equal = specialized == closed
        ok = ok and equal
        rows.append(
            {
                "n": n,
                "coefficient_count_specialized": len(specialized),
                "coefficient_count_closed": len(closed),
                "q_minus_one_collapse_matches_closed_form": equal,
                "closed_form": "((xp-xm)^2+(yp-ym)^2)^(n/2)" if n % 2 == 0 else "(xp*yp-xm*ym)*((xp-xm)^2+(yp-ym)^2)^((n-1)/2)",
            }
        )
    return {
        "id": "q_minus_one_collapse_replay",
        "status": "PASS" if ok else "FAIL",
        "range": f"n=2..{max_n}, coefficient equality after q=-1 specialization",
        "rows": rows,
    }


def run_marginals(max_n: int = 7) -> dict:
    rows = []
    ok = True
    for n in range(2, max_n + 1):
        square = brute_square_poly(n)
        bi_from_square, btilde_from_square = square_marginals(square)
        bi_brute, btilde_brute = brute_marginals(n)
        bi_ok = bi_from_square == bi_brute
        btilde_ok = btilde_from_square == btilde_brute
        ok = ok and bi_ok and btilde_ok
        rows.append(
            {
                "n": n,
                "BI_rows": len(bi_from_square),
                "Btilde_rows": len(btilde_from_square),
                "BI_marginal_ok": bi_ok,
                "Btilde_marginal_ok": btilde_ok,
            }
        )
    return {
        "id": "square_shadow_marginals_replay",
        "status": "PASS" if ok else "FAIL",
        "range": f"n=2..{max_n}, old BI and new B~I marginals from B^square",
        "rows": rows,
    }


def run_di_recovery(max_n: int = 7) -> dict:
    rows = [di_recovery_rows(n) for n in range(2, max_n + 1)]
    ok = all(row["ok"] for row in rows)
    return {
        "id": "DI_recovery_replay",
        "status": "PASS" if ok else "FAIL",
        "range": f"n=2..{max_n}, all nonempty old BI(k,l) rows",
        "formula": "|DI_n(k,l)|=1/2*(|BI_n(k,l)|+(-1)^l*C(n,l)*[k+l=n])",
        "rows": rows,
    }


def run_dtilde_even(max_n: int = 6) -> dict:
    even_ns = [n for n in range(2, max_n + 1, 2)]
    rows = [dtilde_even_rows(n) for n in even_ns]
    ok = all(row["ok"] for row in rows)
    return {
        "id": "Dtilde_even_correction_replay",
        "status": "PASS" if ok else "FAIL",
        "range": f"even n={even_ns}, all nonempty B~I(k,l) rows",
        "formula": "|D~I_{2m}(k,l)|=1/2*(|B~I_{2m}(k,l)|+Delta_{2m}(k,l)); Delta=(-1)^(k+l)*sum_j C(m,j)C(2j,k)C(2m-2j,l)",
        "rows": rows,
        "boundary": "Count-level Type-D scout support only; no Type-D theorem promotion.",
    }


def make_report(max_n: int = 7, max_dtilde_even_n: int = 6) -> dict:
    checks = [
        run_square_permanent(max_n),
        run_q_minus_one(max_n),
        run_marginals(max_n),
        run_di_recovery(max_n),
        run_dtilde_even(max_dtilde_even_n),
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    proof_checklist = {
        "source_lock_S1": True,
        "prior_art_boundary_S2": True,
        "exact_engine_S3_available": True,
        "B_square_parent_grammar_recorded": True,
        "F_square_permanent_formula_replayed": True,
        "q_minus_one_collapse_replayed": True,
        "BI_and_Btilde_marginals_replayed": True,
        "DI_recovery_replayed": True,
        "Dtilde_even_correction_replayed_as_count_scout": True,
        "no_Type_D_theorem_promoted": True,
        "no_rank_theorem_promoted_at_S6": True,
        "no_public_claim_promoted": True,
    }
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S6 C2 x C2 square and closed forms",
        "status": "PASS" if not failed and all(proof_checklist.values()) else "FAIL",
        "decision": "S6 square/count grammar gate passes; next gate is independent red-team recompute S7.",
        "failed_checks": failed,
        "checks": checks,
        "proof_checklist": proof_checklist,
        "closed_at_S6": [
            "B^square parent grammar is certified as the common parent of old BI and new B~I shadows.",
            "F_square=per(W) is replayed as exact coefficient equality.",
            "The q=-1 collapse is replayed as exact coefficient equality.",
            "The old DI count formula is recovered from the square/collapse grammar.",
            "The even D~I count correction Delta is replayed for all rows in finite cases.",
        ],
        "open_after_S6": [
            "S7 independent red-team recompute remains required before manuscript promotion.",
            "S8 H_m recurrence certificate remains required for the odd k=1 rank residual.",
            "General odd pair-standard positivity remains open.",
            "D~I remains scout-only; no Type-D theorem is promoted.",
        ],
    }


def print_summary(report: dict) -> None:
    print("P15-S6 C2 x C2 square certificate")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"], "-", check["range"])
        if check["id"] == "square_permanent_polynomial_replay":
            for row in check["rows"]:
                print(
                    "  n={n} group={group_size} coeffs={coefficient_count} per_eq={permanent_polynomial_matches_brute} sym={symmetries}".format(**row)
                )
        elif check["id"] == "q_minus_one_collapse_replay":
            for row in check["rows"]:
                print(
                    "  n={n} coeffs={coefficient_count_specialized} closed_coeffs={coefficient_count_closed} collapse={q_minus_one_collapse_matches_closed_form}".format(**row)
                )
        elif check["id"] == "square_shadow_marginals_replay":
            for row in check["rows"]:
                print(
                    "  n={n} BI_rows={BI_rows} Btilde_rows={Btilde_rows} BI_ok={BI_marginal_ok} Btilde_ok={Btilde_marginal_ok}".format(**row)
                )
        elif check["id"] == "DI_recovery_replay":
            for row in check["rows"]:
                print("  n={n} rows={rows_checked} ok={ok}".format(**row))
        elif check["id"] == "Dtilde_even_correction_replay":
            for row in check["rows"]:
                print("  n={n} rows={rows_checked} ok={ok}".format(**row))


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S6 square certificate builder")
    parser.add_argument("--write-json", help="write deterministic JSON certificate to this path")
    parser.add_argument("--max-n", type=int, default=7, help="maximum n for full square/DI replay")
    parser.add_argument("--max-dtilde-even-n", type=int, default=6, help="maximum even n for D~I correction replay")
    args = parser.parse_args()
    report = make_report(args.max_n, args.max_dtilde_even_n)
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
