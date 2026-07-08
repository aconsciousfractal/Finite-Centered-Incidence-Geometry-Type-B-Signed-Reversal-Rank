# -*- coding: utf-8 -*-
"""
P15-S11 bounded raw rank smoke check.

This is deliberately independent of the S3/S7 engines. It rebuilds signed
permutations, diagonal B~I_n(k,k) layers, the two channel matrices, and exact
ranks over Q for a small bounded set of rows.

This script is a reviewer-facing smoke check, not a proof engine. It is bounded
at n<=6 to avoid exploratory long runs.

Run with python -B to avoid writing __pycache__ on Windows.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Iterable

SCHEMA = "p15.s11.raw_rank_smoke_check.v1"
RANK_CASES = [(3, 1), (4, 1), (4, 2), (5, 1), (5, 2), (5, 3), (6, 1), (6, 2)]
MAX_BOUNDARY_N = 6


def q(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value, 1)


def zero_matrix(n: int) -> list[list[int]]:
    return [[0 for _ in range(n)] for _ in range(n)]


def rank_exact(matrix: list[list[int | Fraction]]) -> int:
    if not matrix:
        return 0
    mat = [[q(x) for x in row] for row in matrix]
    rows = len(mat)
    cols = len(mat[0])
    r = 0
    for c in range(cols):
        pivot = None
        for i in range(r, rows):
            if mat[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        pv = mat[r][c]
        mat[r] = [x / pv for x in mat[r]]
        for i in range(rows):
            if i != r and mat[i][c] != 0:
                factor = mat[i][c]
                mat[i] = [mat[i][j] - factor * mat[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def matmul(a: list[list[int | Fraction]], b: list[list[int | Fraction]]) -> list[list[Fraction]]:
    rows = len(a)
    mid = len(b)
    cols = len(b[0]) if b else 0
    out = [[Fraction(0, 1) for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = q(a[i][k])
            if aik == 0:
                continue
            for j in range(cols):
                out[i][j] += aik * q(b[k][j])
    return out


def standard_projected_matrix(matrix: list[list[int]]) -> list[list[Fraction]]:
    n = len(matrix)
    projector = [[Fraction(1 if i == j else 0, 1) - Fraction(1, n) for j in range(n)] for i in range(n)]
    return matmul(matmul(projector, matrix), projector)


def reversal(n: int) -> tuple[int, ...]:
    return tuple(n - 1 - i for i in range(n))


def signed_perms(n: int):
    for pi in itertools.permutations(range(n)):
        for eps in itertools.product((1, -1), repeat=n):
            yield pi, eps


def positive_hit_counts(pi: tuple[int, ...], eps: tuple[int, ...], rho: tuple[int, ...]) -> tuple[int, int]:
    k = sum(1 for i in range(len(pi)) if pi[i] == i and eps[i] == 1)
    ell = sum(1 for i in range(len(pi)) if pi[i] == rho[i] and eps[i] == 1)
    return k, ell


def expected_nonempty(n: int, k: int) -> bool:
    if n < 2 or k < 1:
        return False
    if n % 2 == 0:
        m = n // 2
        return (1 <= k < m) or (k == m and m % 2 == 0)
    m = n // 2
    return (1 <= k <= m) or (k == m + 1 and m % 2 == 0)


def init_case(n: int) -> dict:
    return {"size": 0, "m_ref": zero_matrix(n), "m_pair": zero_matrix(n)}


def update_case(case: dict, pi: tuple[int, ...], eps: tuple[int, ...]) -> None:
    case["size"] += 1
    m_ref = case["m_ref"]
    m_pair = case["m_pair"]
    for i, j in enumerate(pi):
        m_ref[j][i] += eps[i]
        m_pair[j][i] += 1


def compute_for_n(n: int, ks: Iterable[int]) -> dict[int, dict]:
    rho = reversal(n)
    cases = {k: init_case(n) for k in sorted(set(ks))}
    for pi, eps in signed_perms(n):
        k, ell = positive_hit_counts(pi, eps, rho)
        if k == ell and k in cases:
            update_case(cases[k], pi, eps)
    return cases


def rank_rows() -> list[dict]:
    by_n: dict[int, list[int]] = defaultdict(list)
    for n, k in RANK_CASES:
        by_n[n].append(k)
    rows = []
    for n in sorted(by_n):
        cases = compute_for_n(n, by_n[n])
        for k in sorted(by_n[n]):
            case = cases[k]
            rank_ref = rank_exact(case["m_ref"])
            rank_pair_std = rank_exact(standard_projected_matrix(case["m_pair"]))
            expected_ref = (n + 1) // 2
            expected_pair_std = expected_ref - 1
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "layer_size": case["size"],
                    "rank_ref": rank_ref,
                    "rank_pair_std": rank_pair_std,
                    "expected_ref": expected_ref,
                    "expected_pair_std": expected_pair_std,
                    "ok": case["size"] > 0 and rank_ref == expected_ref and rank_pair_std == expected_pair_std,
                }
            )
    return rows


def boundary_rows(max_n: int = MAX_BOUNDARY_N) -> list[dict]:
    rows = []
    for n in range(2, max_n + 1):
        top = (n + 2) // 2
        cases = compute_for_n(n, range(1, top + 1))
        for k in range(1, top + 1):
            size = cases[k]["size"]
            expected = expected_nonempty(n, k)
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "layer_size": size,
                    "expected_nonempty": expected,
                    "observed_nonempty": size > 0,
                    "ok": (size > 0) == expected,
                }
            )
    return rows


def make_report() -> dict:
    ranks = rank_rows()
    boundary = boundary_rows()
    checks = [
        {
            "id": "raw_rank_smoke_cases",
            "status": "PASS" if all(row["ok"] for row in ranks) else "FAIL",
            "rows": ranks,
        },
        {
            "id": "nonempty_predicate_small_replay",
            "status": "PASS" if all(row["ok"] for row in boundary) else "FAIL",
            "range": f"n=2..{MAX_BOUNDARY_N}, k up to ceil(n/2)+boundary",
            "rows": boundary,
        },
    ]
    failed = [check["id"] for check in checks if check["status"] != "PASS"]
    return {
        "schema": SCHEMA,
        "project": "P15 signed reversal rank theorem",
        "gate": "P15-S11 bounded raw rank smoke check",
        "status": "PASS" if not failed else "FAIL",
        "decision": "Bounded independent raw rank smoke check passes." if not failed else "Smoke check failed; inspect rows.",
        "failed_checks": failed,
        "scope": "Small exact raw recompute only; not a proof engine and not a public theorem promotion.",
        "rank_cases": RANK_CASES,
        "checks": checks,
    }


def write_markdown(report: dict, path: Path) -> None:
    lines = [
        "# P15-S11 Bounded Raw Rank Smoke Check Certificate",
        "",
        "Date: 2026-07-07",
        f"Status: {report['status']}",
        "",
        "This verifier does not import S3 or S7. It rebuilds signed permutations, diagonal layers, channel matrices, standard projection, and exact ranks over Q for a bounded set of rows.",
        "",
        "## Rank Rows",
        "",
        "```text",
    ]
    for row in report["checks"][0]["rows"]:
        lines.append(
            "n={n} k={k} size={layer_size} ranks=({rank_ref},{rank_pair_std}) expected=({expected_ref},{expected_pair_std}) ok={ok}".format(**row)
        )
    lines.extend([
        "```",
        "",
        "## Nonempty Boundary Replay",
        "",
        "```text",
    ])
    for row in report["checks"][1]["rows"]:
        lines.append(
            "n={n} k={k} size={layer_size} expected_nonempty={expected_nonempty} observed_nonempty={observed_nonempty} ok={ok}".format(**row)
        )
    lines.extend([
        "```",
        "",
        "## Boundary",
        "",
        "This is a compact smoke check only. The theorem still relies on the S4/S5/S8/S9/S9C proof chain and the S11 manuscript lemmas.",
        "",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def print_summary(report: dict) -> None:
    print("P15-S11 bounded raw rank smoke check")
    print("schema:", report["schema"])
    print("status:", report["status"])
    print("decision:", report["decision"])
    for check in report["checks"]:
        print("CHECK", check["id"], check["status"])
        print("  rows:", len(check["rows"]))


def main() -> int:
    parser = argparse.ArgumentParser(description="P15-S11 bounded raw rank smoke check")
    parser.add_argument("--write-json", help="write JSON certificate")
    parser.add_argument("--write-md", help="write Markdown certificate")
    args = parser.parse_args()
    report = make_report()
    print_summary(report)
    if args.write_json:
        path = Path(args.write_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("wrote_json:", str(path))
    if args.write_md:
        write_markdown(report, Path(args.write_md))
        print("wrote_md:", args.write_md)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
