#!/usr/bin/env python3
"""
Exact finite certificate for the native full-B_n zero family

    M_{(1.1 | nu)}(B~I_8(3,3)) = 0,  nu |- 6.

This extends the single zero-block witness without promoting any full-fingerprint
classification theorem.  It reuses the exact little-group machinery from
p15_full_bn_zero_block_certificate.py.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from p15_full_bn_zero_block_certificate import (
    Mask,
    Perm,
    layer_size,
    masks_of_size,
    restricted_perm_data,
    scaled_permutation_matrices,
    sign_of_perm,
    sign_poly_target_coeff,
)


@dataclass(frozen=True)
class ShapeRow:
    lambda_minus: str
    dim: int
    common_denominator: int
    max_abs_scaled_entry: int
    nonzero_scaled_entries: int
    relation_checks_pass: bool
    status: str


@dataclass(frozen=True)
class FamilySummary:
    status: str
    claim: str
    n: int
    k: int
    l: int
    lambda_plus: str
    lambda_minus_count: int
    layer_size: int
    subset_count: int
    active_terms: int
    sum_abs_coefficients: int
    rows: List[ShapeRow]
    runtime_seconds: float
    boundary: str


def partitions(n: int, max_part: int | None = None):
    if n == 0:
        yield ()
        return
    if max_part is None or max_part > n:
        max_part = n
    for first in range(max_part, 0, -1):
        next_max = min(first, n - first) if n - first else 0
        for rest in partitions(n - first, next_max):
            yield (first,) + rest


def part_str(part: Tuple[int, ...]) -> str:
    return "1" if part == (1,) else ".".join(str(x) for x in part)


def precompute_active_terms(n: int, target: Tuple[int, int], plus_size: int):
    subset_masks = masks_of_size(n, plus_size)
    subset_index = {mask: i for i, mask in enumerate(subset_masks)}
    full_mask = (1 << n) - 1
    terms = []
    restrict_cache: Dict[Tuple[Perm, Mask], Tuple[Mask, Perm]] = {}
    sum_abs_coefficients = 0
    for pi in itertools.permutations(range(n)):
        pi = tuple(pi)
        for src_si, A_mask in enumerate(subset_masks):
            B_mask = full_mask ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            data_a = restrict_cache.get((pi, A_mask))
            if data_a is None:
                data_a = restricted_perm_data(pi, A_mask, n)
                restrict_cache[(pi, A_mask)] = data_a
            Ap_mask, sigmaA = data_a
            data_b = restrict_cache.get((pi, B_mask))
            if data_b is None:
                data_b = restricted_perm_data(pi, B_mask, n)
                restrict_cache[(pi, B_mask)] = data_b
            _Bp_mask, sigmaB = data_b
            dst_si = subset_index[Ap_mask]
            coeff_signed = coeff * sign_of_perm(sigmaA)
            terms.append((dst_si, src_si, coeff_signed, sigmaB))
            sum_abs_coefficients += abs(coeff)
    return subset_masks, terms, sum_abs_coefficients


def certify_shape(shape: Tuple[int, ...], subset_count: int, terms) -> ShapeRow:
    denom, minus_mats, checks = scaled_permutation_matrices(shape)
    relation_checks_pass = all(checks.values())
    dim = len(next(iter(minus_mats.values())))
    agg = np.zeros((subset_count, subset_count, dim, dim), dtype=np.int64)
    for dst_si, src_si, coeff_signed, sigmaB in terms:
        agg[dst_si, src_si] += coeff_signed * minus_mats[sigmaB]
    max_abs = int(np.max(np.abs(agg))) if agg.size else 0
    nonzero = int(np.count_nonzero(agg))
    status = "PASS" if relation_checks_pass and max_abs == 0 and nonzero == 0 else "FAIL"
    return ShapeRow(
        lambda_minus=part_str(shape),
        dim=dim,
        common_denominator=denom,
        max_abs_scaled_entry=max_abs,
        nonzero_scaled_entries=nonzero,
        relation_checks_pass=relation_checks_pass,
        status=status,
    )


def write_markdown(summary: FamilySummary, path: Path) -> None:
    lines = [
        "# P15 Full Bn Zero Family Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Method",
        "",
        "The certificate uses the exact little-group model from `p15_full_bn_zero_block_certificate.py`.",
        "The plus shape is fixed at `lambda_plus=(1,1)`, so the plus Specht action is the sign representation of `S_2`.",
        "For each partition `nu` of `6`, the `S^nu` seminormal matrices are computed over `Q`, a common denominator is cleared, and the resulting integer block matrix is summed exactly.",
        "",
        "## Shared Layer Audit",
        "",
        "| field | value |",
        "|---|---:|",
        f"| `layer_size` | {summary.layer_size} |",
        f"| `subset_count` | {summary.subset_count} |",
        f"| `active_terms` | {summary.active_terms} |",
        f"| `sum_abs_coefficients` | {summary.sum_abs_coefficients} |",
        "",
        "## Shape Rows",
        "",
        "| lambda_minus | dim S^nu | denominator | max abs cleared entry | nonzero cleared entries | status |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary.rows:
        lines.append(
            f"| `{row.lambda_minus}` | {row.dim} | {row.common_denominator} | {row.max_abs_scaled_entry} | {row.nonzero_scaled_entries} | {row.status} |"
        )
    lines += [
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "This is a finite exact zero-family witness. It does not promote a full native `B_n` fingerprint theorem or a classification theorem.",
        "",
        "Prossima task: decide whether to seek a conceptual representation-theoretic explanation for this `lambda_plus=(1,1)` zero family, or return to final PDF readthrough/freeze.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run_family(out_json: Path | None, out_md: Path | None) -> FamilySummary:
    t0 = time.time()
    n = 8
    target = (3, 3)
    plus_shape = (1, 1)
    minus_shapes = list(partitions(6))
    subset_masks, terms, sum_abs_coefficients = precompute_active_terms(n, target, sum(plus_shape))
    rows = [certify_shape(shape, len(subset_masks), terms) for shape in minus_shapes]
    status = "PASS" if all(row.status == "PASS" for row in rows) else "FAIL"
    summary = FamilySummary(
        status=status,
        claim="M_{(1.1 | nu)}(B~I_8(3,3)) = 0 for every partition nu of 6, in the rational little-group model",
        n=n,
        k=target[0],
        l=target[1],
        lambda_plus="1.1",
        lambda_minus_count=len(minus_shapes),
        layer_size=layer_size(n, target),
        subset_count=len(subset_masks),
        active_terms=len(terms),
        sum_abs_coefficients=sum_abs_coefficients,
        rows=rows,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite exact full-fingerprint zero-family witness only; not a public full-B_n fingerprint theorem.",
    )
    if out_json is not None:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md is not None:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = run_family(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())