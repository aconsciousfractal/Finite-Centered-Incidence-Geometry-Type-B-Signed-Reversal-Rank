#!/usr/bin/env python3
"""
Pre-Specht fixed-plus decompression scout for the P15 n=8 defect atlas.

Targets:
    T_1^(8,2) in Mat_8(Q[S_7])
    T_2^(8,3) in Mat_28(Q[S_6])

For a fixed one-row plus shape alpha=(a), the plus Specht action is trivial.
The layer operator can therefore be compressed before choosing a complement
Specht module: aggregate coefficients in a matrix over Q[S_{n-a}], then apply
each complement shape nu only afterwards.

Implementation note: the compressed operators have only 64 and 96 complement
permutations respectively, so this script computes only the Specht matrices that
actually occur, not all elements of S_7 or S_6.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import time
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from p15_full_bn_zero_block_certificate import (
    Mask,
    Perm,
    SpechtSeminormalExact,
    layer_size,
    masks_of_size,
    restricted_perm_data,
    sign_poly_target_coeff,
)

ROOT = Path(__file__).resolve().parents[1]
DEFECT_CSV = ROOT / "experiments" / "external_agents_2026_07_07" / "p15_full_bn_fingerprint_n8_kge2_defects.csv"


@dataclass(frozen=True)
class ShapeRow:
    lambda_minus: str
    specht_dim: int
    total_dim: int
    r_plus_bound_mod_p: int
    rank_mod_p: int
    defect_to_bound: int
    bucket: str
    listed_in_defect_csv: bool
    expected_rank_mod_p: int | None
    expected_defect_to_bound: int | None
    expected_bucket: str | None
    csv_match: bool | None


@dataclass(frozen=True)
class CaseSummary:
    label: str
    n: int
    k: int
    lambda_plus: str
    plus_size: int
    complement_size: int
    prime: int
    layer_size: int
    subset_count: int
    permutations_scanned: int
    raw_active_terms: int
    raw_sum_abs_coefficients: int
    aggregate_group_keys: int
    aggregate_zero_keys_removed: int
    distinct_complement_permutations: int
    max_abs_group_coefficient: int
    rows: List[ShapeRow]
    bucket_counts: Dict[str, int]
    defect_csv_mismatches: int
    one_dimensional_exact_checks: Dict[str, Dict[str, int]]


@dataclass(frozen=True)
class ScoutSummary:
    status: str
    claim: str
    cases: List[CaseSummary]
    runtime_seconds: float
    boundary: str


def partitions(n: int, max_part: int | None = None):
    if n == 0:
        yield ()
        return
    if max_part is None or max_part > n:
        max_part = n
    for first in range(max_part, 0, -1):
        rest_n = n - first
        if rest_n == 0:
            yield (first,)
        else:
            for rest in partitions(rest_n, min(first, rest_n)):
                yield (first,) + rest


def part_str(part: Tuple[int, ...]) -> str:
    if not part:
        return "empty"
    return ".".join(str(x) for x in part)


def rho_perm(n: int) -> Perm:
    return tuple(n - 1 - i for i in range(n))


def load_defect_rows() -> Dict[Tuple[int, int, str, str], Dict[str, str]]:
    out: Dict[Tuple[int, int, str, str], Dict[str, str]] = {}
    if not DEFECT_CSV.exists():
        return out
    with DEFECT_CSV.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (int(row["n"]), int(row["k"]), row["lambda_plus"], row["lambda_minus"])
            out[key] = row
    return out


def mod_rank(mat: np.ndarray, p: int) -> int:
    A = np.array(mat % p, dtype=np.int64, copy=True)
    rows, cols = A.shape
    rank = 0
    for col in range(cols):
        nz = np.flatnonzero(A[rank:, col] % p)
        if nz.size == 0:
            continue
        pivot = int(nz[0] + rank)
        if pivot != rank:
            A[[rank, pivot], :] = A[[pivot, rank], :]
        inv = pow(int(A[rank, col] % p), p - 2, p)
        A[rank, :] = (A[rank, :] * inv) % p
        nz_rows = np.flatnonzero(A[:, col] % p)
        for r in nz_rows:
            if int(r) == rank:
                continue
            fac = int(A[r, col] % p)
            if fac:
                A[r, :] = (A[r, :] - fac * A[rank, :]) % p
        rank += 1
        if rank == rows:
            break
    return rank


def fraction_to_mod(x: Fraction, p: int) -> int:
    den = x.denominator % p
    if den == 0:
        raise RuntimeError(f"denominator {x.denominator} is divisible by prime {p}")
    return (x.numerator % p) * pow(den, p - 2, p) % p


def needed_specht_matrices(shape: Tuple[int, ...], perms: set[Perm], p: int) -> Tuple[Dict[Perm, np.ndarray], int]:
    S = SpechtSeminormalExact(shape)
    mats: Dict[Perm, np.ndarray] = {}
    for perm in perms:
        Mq = S.mat_perm(perm)
        Z = np.zeros((S.dim, S.dim), dtype=np.int64)
        for i, row in enumerate(Mq):
            for j, x in enumerate(row):
                if x:
                    Z[i, j] = fraction_to_mod(x, p)
        mats[perm] = Z
    return mats, S.dim


def build_prespecht_terms(n: int, k: int, plus_size: int):
    target = (k, k)
    subset_masks = masks_of_size(n, plus_size)
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full_mask = (1 << n) - 1
    grouped: Dict[Tuple[int, int, Perm], int] = {}
    raw_active = 0
    raw_abs = 0
    permutations_scanned = 0
    zero_keys_removed = 0
    for pi0 in itertools.permutations(range(n)):
        permutations_scanned += 1
        pi = tuple(pi0)
        for src_si, A_mask in enumerate(subset_masks):
            B_mask = full_mask ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            raw_active += 1
            raw_abs += abs(coeff)
            Ap_mask, _sigmaA = restricted_perm_data(pi, A_mask, n)
            _Bp_mask, sigmaB = restricted_perm_data(pi, B_mask, n)
            dst_si = subset_index[Ap_mask]
            key = (dst_si, src_si, sigmaB)
            new_val = grouped.get(key, 0) + coeff
            if new_val:
                grouped[key] = new_val
            elif key in grouped:
                del grouped[key]
                zero_keys_removed += 1
    return subset_masks, grouped, permutations_scanned, raw_active, raw_abs, zero_keys_removed


def r_sigmas(n: int, plus_size: int, subset_masks: List[Mask]) -> set[Perm]:
    full_mask = (1 << n) - 1
    rho = rho_perm(n)
    out: set[Perm] = set()
    for A_mask in subset_masks:
        B_mask = full_mask ^ A_mask
        _Bp_mask, sigmaB = restricted_perm_data(rho, B_mask, n)
        out.add(sigmaB)
    return out


def block_matrix_mod(grouped, subset_count: int, mats_mod: Dict[Perm, np.ndarray], dim: int, p: int) -> np.ndarray:
    M = np.zeros((subset_count * dim, subset_count * dim), dtype=np.int64)
    for (dst_si, src_si, sigmaB), coeff in grouped.items():
        Z = mats_mod[sigmaB]
        r0 = dst_si * dim
        c0 = src_si * dim
        M[r0:r0 + dim, c0:c0 + dim] = (M[r0:r0 + dim, c0:c0 + dim] + coeff * Z) % p
    return M


def r_operator_mod(n: int, plus_size: int, subset_masks: List[Mask], mats_mod: Dict[Perm, np.ndarray], dim: int, p: int) -> np.ndarray:
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full_mask = (1 << n) - 1
    rho = rho_perm(n)
    R = np.zeros((len(subset_masks) * dim, len(subset_masks) * dim), dtype=np.int64)
    for src_si, A_mask in enumerate(subset_masks):
        B_mask = full_mask ^ A_mask
        Ap_mask, _sigmaA = restricted_perm_data(rho, A_mask, n)
        _Bp_mask, sigmaB = restricted_perm_data(rho, B_mask, n)
        dst_si = subset_index[Ap_mask]
        Z = mats_mod[sigmaB]
        r0 = dst_si * dim
        c0 = src_si * dim
        R[r0:r0 + dim, c0:c0 + dim] = (R[r0:r0 + dim, c0:c0 + dim] + Z) % p
    return R


def perm_sign(perm: Perm) -> int:
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            inv += int(perm[i] > perm[j])
    return -1 if inv % 2 else 1


def one_dimensional_exact(grouped, subset_count: int, complement_size: int) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for label, char in [
        (str(complement_size), lambda sigma: 1),
        (".".join(["1"] * complement_size), lambda sigma: perm_sign(sigma)),
    ]:
        M = [[0 for _ in range(subset_count)] for __ in range(subset_count)]
        for (dst, src, sigma), coeff in grouped.items():
            M[dst][src] += coeff * char(sigma)
        vals = [x for row in M for x in row]
        out[label] = {
            "nonzero_entries": sum(1 for x in vals if x),
            "max_abs_entry": max([abs(x) for x in vals] or [0]),
            "rank_mod_1000003": mod_rank(np.array(M, dtype=np.int64), 1000003),
        }
    return out


def run_case(n: int, k: int, plus_size: int, prime: int, defect_rows) -> CaseSummary:
    lambda_plus = str(plus_size)
    subset_masks, grouped, permutations_scanned, raw_active, raw_abs, zero_removed = build_prespecht_terms(n, k, plus_size)
    complement_size = n - plus_size
    needed = {sigma for (_dst, _src, sigma) in grouped.keys()} | r_sigmas(n, plus_size, subset_masks)
    rows: List[ShapeRow] = []
    bucket_counts: Dict[str, int] = {}
    mismatches = 0
    I_cache: Dict[int, np.ndarray] = {}
    for shape in partitions(complement_size):
        lambda_minus = part_str(shape)
        mats_mod, dim = needed_specht_matrices(shape, needed, prime)
        M = block_matrix_mod(grouped, len(subset_masks), mats_mod, dim, prime)
        rank = mod_rank(M, prime)
        R = r_operator_mod(n, plus_size, subset_masks, mats_mod, dim, prime)
        total_dim = len(subset_masks) * dim
        I = I_cache.get(total_dim)
        if I is None:
            I = np.eye(total_dim, dtype=np.int64)
            I_cache[total_dim] = I
        bound = mod_rank((I + R) % prime, prime)
        defect = bound - rank
        bucket = "ZERO" if rank == 0 and bound > 0 else ("EXTRA_DEFECT" if defect > 0 else "EQUAL_R_BOUND")
        bucket_counts[bucket] = bucket_counts.get(bucket, 0) + 1
        expected = defect_rows.get((n, k, lambda_plus, lambda_minus))
        listed = expected is not None
        exp_rank = int(expected["rank_mod_p"]) if expected else None
        exp_defect = int(expected["defect_to_bound"]) if expected else None
        exp_bucket = expected["bucket"] if expected else None
        csv_match = None
        if expected:
            csv_match = (exp_rank == rank and exp_defect == defect and exp_bucket == bucket)
            if not csv_match:
                mismatches += 1
        elif bucket != "EQUAL_R_BOUND":
            csv_match = False
            mismatches += 1
        rows.append(ShapeRow(
            lambda_minus=lambda_minus,
            specht_dim=dim,
            total_dim=total_dim,
            r_plus_bound_mod_p=bound,
            rank_mod_p=rank,
            defect_to_bound=defect,
            bucket=bucket,
            listed_in_defect_csv=listed,
            expected_rank_mod_p=exp_rank,
            expected_defect_to_bound=exp_defect,
            expected_bucket=exp_bucket,
            csv_match=csv_match,
        ))
    return CaseSummary(
        label=f"T_{lambda_plus}^(8,{k})",
        n=n,
        k=k,
        lambda_plus=lambda_plus,
        plus_size=plus_size,
        complement_size=complement_size,
        prime=prime,
        layer_size=layer_size(n, (k, k)),
        subset_count=len(subset_masks),
        permutations_scanned=permutations_scanned,
        raw_active_terms=raw_active,
        raw_sum_abs_coefficients=raw_abs,
        aggregate_group_keys=len(grouped),
        aggregate_zero_keys_removed=zero_removed,
        distinct_complement_permutations=len(needed),
        max_abs_group_coefficient=max([abs(x) for x in grouped.values()] or [0]),
        rows=rows,
        bucket_counts=dict(sorted(bucket_counts.items())),
        defect_csv_mismatches=mismatches,
        one_dimensional_exact_checks=one_dimensional_exact(grouped, len(subset_masks), complement_size),
    )


def write_markdown(summary: ScoutSummary, path: Path) -> None:
    lines = [
        "# P15 Pre-Specht Fixed-Plus Decompression Scout",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Purpose",
        "",
        "This scout fixes a small plus shape before applying complement Specht modules. It tests whether whole families of n=8 full-fingerprint defects are better explained by a compressed operator over a complement group algebra.",
        "",
        "```text",
        "T_1^(8,2) in Mat_8(Q[S_7])",
        "T_2^(8,3) in Mat_28(Q[S_6])",
        "```",
        "",
        "The ranks below are modular over the listed prime, with the `R=+` bound recomputed from the same little-group model. Rows are compared against the external n=8 defect CSV when the CSV lists them.",
        "",
    ]
    for case in summary.cases:
        lines += [
            f"## {case.label}",
            "",
            "| field | value |",
            "|---|---:|",
            f"| `layer_size` | {case.layer_size} |",
            f"| `subset_count` | {case.subset_count} |",
            f"| `raw_active_terms` | {case.raw_active_terms} |",
            f"| `aggregate_group_keys` | {case.aggregate_group_keys} |",
            f"| `distinct_complement_permutations` | {case.distinct_complement_permutations} |",
            f"| `max_abs_group_coefficient` | {case.max_abs_group_coefficient} |",
            f"| `defect_csv_mismatches` | {case.defect_csv_mismatches} |",
            "",
            "Bucket counts:",
            "",
            "```text",
            json.dumps(case.bucket_counts, sort_keys=True),
            "```",
            "",
            "One-dimensional complement checks:",
            "",
            "```text",
            json.dumps(case.one_dimensional_exact_checks, indent=2, sort_keys=True),
            "```",
            "",
            "| lambda_minus | dim S^nu | total dim | R+ bound | rank mod p | defect | bucket | CSV match |",
            "|---|---:|---:|---:|---:|---:|---|---|",
        ]
        for row in case.rows:
            match = "NA" if row.csv_match is None else str(row.csv_match)
            lines.append(
                f"| `{row.lambda_minus}` | {row.specht_dim} | {row.total_dim} | {row.r_plus_bound_mod_p} | {row.rank_mod_p} | {row.defect_to_bound} | {row.bucket} | {match} |"
            )
        lines.append("")
    lines += [
        "## Reading",
        "",
        "`T_1^(8,2)` accounts for a full fixed-plus family: every complement partition of 7 is defective, with two exact zero quotients at the trivial and sign complement characters. Downstream certificates now explain the two zero rows, the natural quotient layers, the sign-twist conjugates, and the three residual Specht quotient mechanisms; the finite `lambda_plus=1`, `n=8,k=2` family is fully accounted.",
        "",
        "`T_2^(8,3)` accounts for a second compact family: ten of the eleven complement partitions of 6 are extra-defective, with no zero row in this fixed-plus family. This remains the next fixed-plus extra-defect target after the `T_1^(8,2)` family closure.",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: keep this package appendix-scoped; optional rational audits and any remaining families are separate finite checks, not a classification theorem.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run(out_json: Path | None, out_md: Path | None) -> ScoutSummary:
    t0 = time.time()
    defect_rows = load_defect_rows()
    prime = 1000003
    cases = [
        run_case(8, 2, 1, prime, defect_rows),
        run_case(8, 3, 2, prime, defect_rows),
    ]
    status = "PASS" if all(case.defect_csv_mismatches == 0 for case in cases) else "FAIL"
    summary = ScoutSummary(
        status=status,
        claim="Pre-Specht fixed-plus decompression reproduces the n=8 defect CSV families for lambda_plus=1 at k=2 and lambda_plus=2 at k=3.",
        cases=cases,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite modular scout and family-decompression artifact only; no full native B_n fingerprint theorem or classification theorem.",
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = run(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
