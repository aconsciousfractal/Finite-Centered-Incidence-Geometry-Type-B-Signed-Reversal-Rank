#!/usr/bin/env python3
"""
Exact rational direct-rank audit for one-row fixed-plus packages T_1, T_2, T_4.

The modular quotient certificates use ranks over p=1000003.  This audit recomputes
the direct Specht block rank and the R+ bound over Q for one selected one-row
fixed-plus family at a time.  Seminormal Specht matrices are computed exactly,
scaled to integer matrices, and ranked with SymPy over the rationals.

Run family by family; do not batch all three unless runtime is acceptable:

    --family T4   means lambda_plus=4, n=8,k=3
    --family T2   means lambda_plus=2, n=8,k=3
    --family T1   means lambda_plus=1, n=8,k=2
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from sympy.polys.domains import ZZ
from sympy.polys.matrices import DomainMatrix

from p15_full_bn_zero_block_certificate import SpechtSeminormalExact, restricted_perm_data, rho_perm
from p15_full_bn_prespecht_fixed_plus_scout import (
    build_prespecht_terms,
    load_defect_rows,
    part_str,
    partitions,
    r_sigmas,
)

ROOT = Path(__file__).resolve().parents[1]


def gcd_lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


@dataclass(frozen=True)
class RationalRow:
    lambda_minus: str
    specht_dim: int
    total_dim: int
    common_denominator: int
    rank_Q: int
    r_plus_bound_Q: int
    defect_to_bound_Q: int
    bucket_Q: str
    listed_in_defect_csv: bool
    expected_rank_mod_p: int | None
    expected_defect_to_bound: int | None
    expected_bucket: str | None
    csv_match: bool | None


@dataclass(frozen=True)
class RationalFamilySummary:
    status: str
    claim: str
    family: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    subset_count: int
    group_keys: int
    raw_active_terms: int
    raw_sum_abs_coefficients: int
    zero_keys_removed: int
    rows: List[RationalRow]
    runtime_seconds: float
    boundary: str
    next_task: str


def needed_scaled_specht_matrices(shape: Tuple[int, ...], perms: set[Tuple[int, ...]]):
    S = SpechtSeminormalExact(shape)
    mats_q = {}
    denom = 1
    for perm in perms:
        M = S.mat_perm(perm)
        mats_q[perm] = M
        for row in M:
            for x in row:
                denom = gcd_lcm(denom, x.denominator)
    mats: Dict[Tuple[int, ...], np.ndarray] = {}
    for perm, M in mats_q.items():
        Z = np.zeros((S.dim, S.dim), dtype=object)
        for i, row in enumerate(M):
            for j, x in enumerate(row):
                y = x * denom
                if y.denominator != 1:
                    raise AssertionError((shape, perm, i, j, x, denom))
                Z[i, j] = int(y.numerator)
        mats[perm] = Z
    return denom, mats, S.dim


def exact_rank_int_matrix(A: np.ndarray) -> int:
    rows = [[int(x) for x in row] for row in A.tolist()]
    return int(DomainMatrix.from_list(rows, ZZ).rank())


def block_matrix_scaled(grouped, subset_count: int, mats: Dict[Tuple[int, ...], np.ndarray], dim: int) -> np.ndarray:
    M = np.zeros((subset_count * dim, subset_count * dim), dtype=object)
    for (dst, src, sigma), coeff in grouped.items():
        r0 = dst * dim
        c0 = src * dim
        M[r0:r0 + dim, c0:c0 + dim] += int(coeff) * mats[sigma]
    return M


def r_operator_scaled(n: int, plus_size: int, plus_masks: List[int], mats: Dict[Tuple[int, ...], np.ndarray], dim: int) -> np.ndarray:
    subset_index = {mask: i for i, mask in enumerate(plus_masks)}
    full_mask = (1 << n) - 1
    rho = rho_perm(n)
    R = np.zeros((len(plus_masks) * dim, len(plus_masks) * dim), dtype=object)
    for src, A_mask in enumerate(plus_masks):
        B_mask = full_mask ^ A_mask
        Ap_mask, _sigma_a = restricted_perm_data(rho, A_mask, n)
        _Bp_mask, sigma_b = restricted_perm_data(rho, B_mask, n)
        dst = subset_index[Ap_mask]
        r0 = dst * dim
        c0 = src * dim
        R[r0:r0 + dim, c0:c0 + dim] += mats[sigma_b]
    return R


def family_params(family: str) -> Tuple[int, int, int, str]:
    if family == "T1":
        return 8, 2, 1, "T_1^(8,2)"
    if family == "T2":
        return 8, 3, 2, "T_2^(8,3)"
    if family == "T4":
        return 8, 3, 4, "T_4^(8,3)"
    raise ValueError(f"unknown family {family!r}")


def parse_shape_filter(value: str | None) -> set[str] | None:
    if not value:
        return None
    return {item.strip() for item in value.split(',') if item.strip()}


def certify_family(family: str, out_json: Path | None, out_md: Path | None, only_shapes: set[str] | None = None) -> RationalFamilySummary:
    t0 = time.time()
    n, k, plus_size, label = family_params(family)
    plus_masks, grouped, _perms, raw_active, raw_abs, zero_removed = build_prespecht_terms(n, k, plus_size)
    needed = {sigma for (_dst, _src, sigma) in grouped} | r_sigmas(n, plus_size, plus_masks)
    defect_rows = load_defect_rows()
    rows: List[RationalRow] = []
    for shape in partitions(n - plus_size):
        lambda_minus = part_str(shape)
        if only_shapes is not None and lambda_minus not in only_shapes:
            continue
        denom, mats, dim = needed_scaled_specht_matrices(shape, needed)
        M = block_matrix_scaled(grouped, len(plus_masks), mats, dim)
        rank_q = exact_rank_int_matrix(M)
        R = r_operator_scaled(n, plus_size, plus_masks, mats, dim)
        bound_q = exact_rank_int_matrix(denom * np.eye(M.shape[0], dtype=object) + R)
        defect = bound_q - rank_q
        bucket = "ZERO" if rank_q == 0 and bound_q > 0 else ("EXTRA_DEFECT" if defect > 0 else "EQUAL_R_BOUND")
        expected = defect_rows.get((n, k, str(plus_size), lambda_minus))
        listed = expected is not None
        exp_rank = int(expected["rank_mod_p"]) if expected else None
        exp_defect = int(expected["defect_to_bound"]) if expected else None
        exp_bucket = expected["bucket"] if expected else None
        csv_match = None if expected is None else (rank_q == exp_rank and defect == exp_defect and bucket == exp_bucket)
        rows.append(RationalRow(
            lambda_minus=lambda_minus,
            specht_dim=dim,
            total_dim=M.shape[0],
            common_denominator=denom,
            rank_Q=rank_q,
            r_plus_bound_Q=bound_q,
            defect_to_bound_Q=defect,
            bucket_Q=bucket,
            listed_in_defect_csv=listed,
            expected_rank_mod_p=exp_rank,
            expected_defect_to_bound=exp_defect,
            expected_bucket=exp_bucket,
            csv_match=csv_match,
        ))
    status = "PASS" if all(row.csv_match is not False for row in rows) else "FAIL"
    summary = RationalFamilySummary(
        status=status,
        claim=f"{label} direct Specht block ranks and R+ bounds have been recomputed exactly over Q; all selected listed defect CSV rows match the rational ranks.",
        family=family,
        n=n,
        k=k,
        lambda_plus=str(plus_size),
        complement_size=n - plus_size,
        subset_count=len(plus_masks),
        group_keys=len(grouped),
        raw_active_terms=raw_active,
        raw_sum_abs_coefficients=raw_abs,
        zero_keys_removed=zero_removed,
        rows=rows,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite one-row fixed-plus rational direct-rank audit only; quotient-layer mechanism certificates remain separate, and no all-n full B_n fingerprint theorem is claimed.",
        next_task="Run the same rational direct-rank audit on the next family, in order T_4, T_2, then T_1; after all pass, update the P15 ledgers from modular to rational-direct where appropriate.",
    )
    if out_json is not None:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md is not None:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: RationalFamilySummary, path: Path) -> None:
    lines = [
        f"# P15 {summary.family} Rational Direct-Rank Audit",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Rows",
        "",
        "| lambda_minus | total dim | denom | rank_Q | R+ bound_Q | defect_Q | bucket_Q | listed CSV | CSV match |",
        "|---|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for row in summary.rows:
        lines.append(
            f"| `{row.lambda_minus}` | {row.total_dim} | {row.common_denominator} | {row.rank_Q} | "
            f"{row.r_plus_bound_Q} | {row.defect_to_bound_Q} | {row.bucket_Q} | {row.listed_in_defect_csv} | {row.csv_match} |"
        )
    lines += [
        "",
        "## Method",
        "",
        "The audit computes Young-seminormal Specht matrices over `Q`, takes a common denominator for the permutations that occur in the compressed fixed-plus operator and in the reversal operator, scales to integer matrices, and calls exact SymPy rank. For `R+`, the matrix ranked is `d I + R_scaled`, which has the same rank as `I + R` over `Q`.",
        "",
        "## Sources",
        "",
        "```text",
        "scripts/p15_full_bn_t124_rational_direct_rank_audit.py",
        "experiments/external_agents_2026_07_07/p15_full_bn_fingerprint_n8_kge2_defects.csv",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        f"Prossima task: {summary.next_task}",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=["T1", "T2", "T4"], required=True)
    parser.add_argument("--write-json", default=None)
    parser.add_argument("--write-md", default=None)
    parser.add_argument("--shapes", default=None, help="Comma-separated lambda_minus shapes, e.g. 7,6.1,5.2")
    args = parser.parse_args()
    summary = certify_family(args.family, Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None, parse_shape_filter(args.shapes))
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
