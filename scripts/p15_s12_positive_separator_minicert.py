#!/usr/bin/env python3
"""S12 mini-certificate for the P15 positive-only separator hunt.

This script is the finite minicertificate behind the promoted S12 Gate A/B
proof notes. It provides an exact rook/inclusion-exclusion count formula for
X_n = B~I_n(1,1), finite cross-checks against the DP scout, finite signed-poset
witness coverage, and exact quotient-window evidence for the translate-tiler
obstruction.
"""
from __future__ import annotations

import argparse
import json
from functools import lru_cache
from math import factorial, gcd
from pathlib import Path

import p15_positive_only_separator_hunt as scout

SCHEMA = "p15.s12_positive_separator_minicert.v1"
EVEN_BLOCK_TERMS = (
    (0, 0, 0, 1),  # empty
    (1, 0, 1, 2),  # one fixed edge in a 2x2 block
    (0, 1, 1, 2),  # one reversal edge in a 2x2 block
    (2, 0, 2, 1),  # both fixed edges in a 2x2 block
    (0, 2, 2, 1),  # both reversal edges in a 2x2 block
)
CENTER_TERMS = (
    (0, 0, 0, 1),  # center not selected
    (1, 0, 1, 1),  # center selected as fixed event
    (0, 1, 1, 1),  # center selected as reversal event
    (1, 1, 1, 1),  # same positive center selected as both events
)


def multiply_terms(poly: dict[tuple[int, int, int], int], terms: tuple[tuple[int, int, int, int], ...]) -> dict[tuple[int, int, int], int]:
    out: dict[tuple[int, int, int], int] = {}
    for (p, q, u), coeff in poly.items():
        for dp, dq, du, weight in terms:
            key = (p + dp, q + dq, u + du)
            out[key] = out.get(key, 0) + coeff * weight
    return out


@lru_cache(maxsize=None)
def rook_poly_coeffs(n: int) -> tuple[tuple[int, int, int, int], ...]:
    poly: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    for _ in range(n // 2):
        poly = multiply_terms(poly, EVEN_BLOCK_TERMS)
    if n % 2:
        poly = multiply_terms(poly, CENTER_TERMS)
    return tuple((p, q, u, coeff) for (p, q, u), coeff in sorted(poly.items()))


@lru_cache(maxsize=None)
def count_x_rook(n: int) -> int:
    total = 0
    for p, q, u, coeff in rook_poly_coeffs(n):
        if p == 0 or q == 0:
            continue
        sign = 1 if (p + q) % 2 == 0 else -1
        total += sign * p * q * coeff * (2 ** (n - u)) * factorial(n - u)
    return total


def quotient_window_target(n: int) -> tuple[int, int, str]:
    if n == 4:
        return 9, 10, "finite_n4"
    if n % 2:
        return 10, 11, "odd_n_ge_5"
    if n <= 20:
        return 11, 12, "even_finite_6_to_20"
    return 10, 11, "even_n_ge_22_target"


def quotient_row(n: int) -> dict:
    x_size = count_x_rook(n)
    group_size = (2**n) * factorial(n)
    lower, upper, label = quotient_window_target(n)
    remainder = group_size % x_size
    return {
        "n": n,
        "X_size": x_size,
        "B_n_size": group_size,
        "quotient_floor": group_size // x_size,
        "B_n_mod_X": remainder,
        "gcd_B_n_X": gcd(group_size, x_size),
        "divides_B_n": remainder == 0,
        "window_label": label,
        "window_lower": lower,
        "window_upper": upper,
        "window_ok": lower * x_size < group_size < upper * x_size,
    }


def rook_vs_dp_rows(n_min: int, n_max: int) -> list[dict]:
    rows = []
    for n in range(n_min, n_max + 1):
        rook = count_x_rook(n)
        dp = scout.count_x_dp(n)
        rows.append({"n": n, "rook_count": rook, "dp_count": dp, "match": rook == dp})
    return rows


def witness_coverage_row(n: int) -> dict:
    labels = tuple(range(-n, 0)) + tuple(range(1, n + 1))
    target = {(a, b) for a in labels for b in labels if a != b}
    covered: set[tuple[int, int]] = set()
    x_size = 0
    for w in scout.signed_permutations(n):
        if not scout.in_x(w):
            continue
        x_size += 1
        order = tuple(v for v in scout.signed_total_order(w) if v != 0)
        for i, earlier in enumerate(order):
            for later in order[i + 1 :]:
                covered.add((later, earlier))
    missing = sorted(target - covered)
    return {
        "n": n,
        "X_size": x_size,
        "ordered_nonzero_pairs": len(target),
        "covered_reversed_pairs": len(covered & target),
        "missing_reversed_pairs": len(missing),
        "all_nonzero_pairs_reversible": not missing,
        "missing_sample": missing[:12],
    }


def make_payload(max_window_n: int, table_window_n: int) -> dict:
    dp_rows = rook_vs_dp_rows(3, 12)
    witness_rows = [witness_coverage_row(n) for n in range(4, 8)]
    quotient_rows_full = [quotient_row(n) for n in range(4, max_window_n + 1)]
    quotient_rows_table = [r for r in quotient_rows_full if r["n"] <= table_window_n or r["n"] in (20, 21, 22)]
    n22 = next((r for r in quotient_rows_full if r["n"] == 22), None)
    return {
        "schema": SCHEMA,
        "object": "X_n = B~I_n(1,1) subset B_n",
        "status": "RESEARCH_MINICERT_NOT_PROMOTED",
        "rook_formula": {
            "even": "R_{2m}(A,B,Z)=(1+2AZ+2BZ+A^2Z^2+B^2Z^2)^m",
            "odd": "R_{2m+1}(A,B,Z)=R_{2m}(A,B,Z)(1+AZ+BZ+ABZ)",
            "count": "N_n=sum_{p,q>=1,u} (-1)^(p+q) p q [A^p B^q Z^u] R_n(A,B,Z) 2^(n-u)(n-u)!",
        },
        "rook_vs_dp_rows_n3_n12": dp_rows,
        "signed_poset_witness_rows_n4_n7": witness_rows,
        "quotient_rows_n4_to_max": quotient_rows_full,
        "quotient_rows_table": quotient_rows_table,
        "summary": {
            "max_window_n": max_window_n,
            "all_rook_dp_matches_n3_n12": all(r["match"] for r in dp_rows),
            "all_nonzero_pair_witness_coverage_n4_n7": all(r["all_nonzero_pairs_reversible"] for r in witness_rows),
            "all_divisibility_fails_n4_to_max": all(not r["divides_B_n"] for r in quotient_rows_full),
            "all_assigned_windows_ok_n4_to_max": all(r["window_ok"] for r in quotient_rows_full),
            "n22_even_switch_floor": n22["quotient_floor"] if n22 else None,
            "n22_even_switch_window_ok": n22["window_ok"] if n22 else None,
            "promotion_gate": "CLOSED_BY_GATE_A_AND_GATE_B: finite minicert retained as support artifact",
        },
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8", newline="\n")


def make_markdown(payload: dict) -> str:
    lines = [
        "# P15 S12 Positive-Only Separator Minicertificate",
        "",
        f"Schema: `{payload['schema']}`",
        "Status: **SUPPORT-MINICERT**, superseded by closed Gate A and Gate B proof notes.",
        "",
        "## Rook Count Formula",
        "",
        f"Even case: `{payload['rook_formula']['even']}`.",
        f"Odd case: `{payload['rook_formula']['odd']}`.",
        f"Count extraction: `{payload['rook_formula']['count']}`.",
        "",
        "Here `A` marks selected fixed-positive events, `B` marks selected reversal-positive events, and `Z` marks distinct forced cells. The `ABZ` center term is the odd-n case where one positive center cell contributes both events.",
        "",
        "## Rook Formula vs DP",
        "",
        "| n | rook count | DP count | match? |",
        "|---:|---:|---:|---|",
    ]
    for r in payload["rook_vs_dp_rows_n3_n12"]:
        lines.append(f"| {r['n']} | {r['rook_count']} | {r['dp_count']} | {r['match']} |")
    lines += [
        "",
        "## Signed-Poset Witness Coverage",
        "",
        "For each ordered nonzero signed pair `(a,b)`, the check asks whether some `w in X_n` puts `b` before `a` in `-w_n<...<-w_1<0<w_1<...<w_n`.",
        "",
        "| n | |X_n| | ordered pairs | covered | missing | all covered? |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for r in payload["signed_poset_witness_rows_n4_n7"]:
        lines.append(
            f"| {r['n']} | {r['X_size']} | {r['ordered_nonzero_pairs']} | "
            f"{r['covered_reversed_pairs']} | {r['missing_reversed_pairs']} | {r['all_nonzero_pairs_reversible']} |"
        )
    lines += [
        "",
        "## Quotient Windows",
        "",
        "The table records exact integer data for `G_n=|B_n|=2^n n!` and `N_n=|X_n|`. The window column is a proof target; this script verifies it exactly through the configured finite range.",
        "",
        "| n | N_n | floor(G_n/N_n) | G_n mod N_n | window | ok? | divides? |",
        "|---:|---:|---:|---:|---|---|---|",
    ]
    for r in payload["quotient_rows_table"]:
        window = f"{r['window_lower']} < G/N < {r['window_upper']} ({r['window_label']})"
        lines.append(
            f"| {r['n']} | {r['X_size']} | {r['quotient_floor']} | {r['B_n_mod_X']} | "
            f"{window} | {r['window_ok']} | {r['divides_B_n']} |"
        )
    lines += [
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Interpretation",
        "",
        "The signed-poset side is closed by the Gate A sign-extension proof and the exact `n=4` witness certificate, including relations with zero.",
        "",
        "The tiler side is closed by the Gate B quotient-window proof. The even quotient window changes at `n=22`, so the earlier `11<G/N<12` heuristic must not be promoted as an all-n statement.",
        "",
        "Promotion status: S12 is now a P11-style separator theorem for `n>=4`; `n=3` is recorded separately as the exact-tiling exception.",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-window-n", type=int, default=80)
    parser.add_argument("--table-window-n", type=int, default=30)
    parser.add_argument("--write-json", type=Path)
    parser.add_argument("--write-md", type=Path)
    args = parser.parse_args()
    payload = make_payload(args.max_window_n, args.table_window_n)
    if args.write_json:
        write_text(args.write_json, json.dumps(payload, indent=2))
    if args.write_md:
        write_text(args.write_md, make_markdown(payload))
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()