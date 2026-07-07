#!/usr/bin/env python3
"""
Scout the lambda_plus=(1,1) group-algebra cancellation for diagonal B~I_n(k,k).

This is deliberately weaker and faster than a full native B_n fingerprint.  It
works before applying any complement Specht module: if the signed layer operator
is zero in the complement group algebra Q[S_{n-2}], then every block
M_{(1.1 | nu)} is zero for nu |- n-2.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from p15_full_bn_zero_block_certificate import (
    Mask,
    Perm,
    layer_size,
    masks_of_size,
    positions_of_mask,
    restricted_perm_data,
    sign_of_perm,
    sign_poly_target_coeff,
)


@dataclass(frozen=True)
class Row:
    n: int
    k: int
    layer_size: int
    subset_count: int
    permutations_scanned: int
    active_terms: int
    group_algebra_keys: int
    nonzero_group_algebra_keys: int
    cancelling_pairs: int
    bad_pairs: int
    total_abs_before_cancellation: int
    total_abs_after_cancellation: int
    abs_unsigned_coefficients: List[int]
    signed_coefficient_patterns: Dict[str, int]
    status: str


@dataclass(frozen=True)
class ScoutSummary:
    status: str
    claim: str
    max_n: int
    rows: List[Row]
    runtime_seconds: float
    boundary: str


def rho(n: int, i: int) -> int:
    return n - 1 - i


def certify_row(n: int, k: int) -> Row:
    target = (k, k)
    plus_size = 2
    subset_masks = masks_of_size(n, plus_size)
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full = (1 << n) - 1
    grouped = defaultdict(list)
    permutations_scanned = 0

    for pi0 in itertools.permutations(range(n)):
        permutations_scanned += 1
        pi = tuple(pi0)
        for src_si, A_mask in enumerate(subset_masks):
            B_mask = full ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            Ap_mask, sigmaA = restricted_perm_data(pi, A_mask, n)
            _Bp_mask, sigmaB = restricted_perm_data(pi, B_mask, n)
            dst_si = subset_index[Ap_mask]
            signed_coeff = coeff * sign_of_perm(sigmaA)
            grouped[(src_si, dst_si, sigmaB)].append((signed_coeff, coeff, sign_of_perm(sigmaA), A_mask, pi))

    nonzero_keys = 0
    cancelling_pairs = 0
    bad_pairs = 0
    total_abs_before = 0
    total_abs_after = 0
    abs_coeffs = set()
    signed_patterns = Counter()

    for (_src_si, _dst_si, _sigmaB), vals in grouped.items():
        signed_sum = sum(v[0] for v in vals)
        total_abs_after += abs(signed_sum)
        total_abs_before += sum(abs(v[0]) for v in vals)
        abs_coeffs.update(abs(v[1]) for v in vals)
        signed_patterns[tuple(sorted(v[0] for v in vals))] += 1
        if signed_sum != 0:
            nonzero_keys += 1
        if len(vals) == 2 and vals[0][0] + vals[1][0] == 0:
            v0, v1 = vals
            A_mask = v0[3]
            A_pts = positions_of_mask(n, A_mask)
            comp_pts = [i for i in range(n) if not ((A_mask >> i) & 1)]
            pi0 = v0[4]
            pi1 = v1[4]
            same_complement = all(pi0[i] == pi1[i] for i in comp_pts)
            swapped_plus = pi0[A_pts[0]] == pi1[A_pts[1]] and pi0[A_pts[1]] == pi1[A_pts[0]]
            same_unsigned = v0[1] == v1[1]
            opposite_plus_sign = v0[2] == -v1[2]
            if same_complement and swapped_plus and same_unsigned and opposite_plus_sign:
                cancelling_pairs += 1
            else:
                bad_pairs += 1
        else:
            bad_pairs += 1

    status = "ZERO_PAIR_CANCEL" if grouped and nonzero_keys == 0 and bad_pairs == 0 else "NONZERO_OR_MIXED"
    return Row(
        n=n,
        k=k,
        layer_size=layer_size(n, target),
        subset_count=len(subset_masks),
        permutations_scanned=permutations_scanned,
        active_terms=sum(len(v) for v in grouped.values()),
        group_algebra_keys=len(grouped),
        nonzero_group_algebra_keys=nonzero_keys,
        cancelling_pairs=cancelling_pairs,
        bad_pairs=bad_pairs,
        total_abs_before_cancellation=total_abs_before,
        total_abs_after_cancellation=total_abs_after,
        abs_unsigned_coefficients=sorted(abs_coeffs),
        signed_coefficient_patterns={str(key): value for key, value in sorted(signed_patterns.items(), key=lambda kv: str(kv[0]))},
        status=status,
    )


def run_scout(max_n: int, out_json: Path | None, out_md: Path | None) -> ScoutSummary:
    t0 = time.time()
    rows: List[Row] = []
    for n in range(4, max_n + 1):
        for k in range(1, n // 2 + 1):
            rows.append(certify_row(n, k))
    summary = ScoutSummary(
        status="PASS",
        claim="Scout for lambda_plus=(1,1) cancellation in Q[S_{n-2}] for diagonal B~I_n(k,k).",
        max_n=max_n,
        rows=rows,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite scout only; it detects exact group-algebra zero rows but does not prove an all-n classification.",
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: ScoutSummary, path: Path) -> None:
    lines = [
        "# P15 lambda_plus=(1,1) Group-Algebra Scout",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Purpose",
        "",
        "This scout tests whether the diagonal layer operator already vanishes in the complement group algebra `Q[S_{n-2}]` for `lambda_plus=(1,1)`. A zero row implies every native block `(1.1|nu)`, `nu |- n-2`, is zero.",
        "",
        "## Rows",
        "",
        "| n | k | layer size | active terms | keys | nonzero keys | cancelling pairs | bad pairs | abs coeffs | status |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in summary.rows:
        lines.append(
            f"| {row.n} | {row.k} | {row.layer_size} | {row.active_terms} | {row.group_algebra_keys} | {row.nonzero_group_algebra_keys} | {row.cancelling_pairs} | {row.bad_pairs} | `{row.abs_unsigned_coefficients}` | {row.status} |"
        )
    zero_rows = [(row.n, row.k) for row in summary.rows if row.status == "ZERO_PAIR_CANCEL"]
    lines += [
        "",
        "## Reading",
        "",
        f"Zero rows detected: `{zero_rows}`.",
        "",
        "The n=8, k=3 row is the only zero-pair-cancellation row in this finite range if the list above is `[(8, 3)]`. In that case the certified zero family is a real special feature of the n=8 middle-minus-one layer, not a generic `lambda_plus=(1,1)` vanishing.",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: either try to prove the finite n=8 pair-cancellation lemma in manuscript style, or extend the scout one notch only if there is a concrete pattern to test.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=8)
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = run_scout(args.max_n, Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
