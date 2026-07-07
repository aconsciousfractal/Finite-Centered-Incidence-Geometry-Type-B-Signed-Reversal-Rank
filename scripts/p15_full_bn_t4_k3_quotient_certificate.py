#!/usr/bin/env python3
"""
Quotient-layer certificate for T_4^(8,3).

This handles the one-row fixed-plus family lambda_plus=4 in the n=8,k=3
full native B_n fingerprint atlas.  The compressed operator is

    T_4^(8,3) in Mat_70(Q[S_4]).

All extra-defect complement rows are accounted for by natural complement-fiber
quotient layers. The remaining sign complement row is equal to its R=+ bound.

This is a finite n=8,k=3 certificate. It is not an all-n theorem and not a
full native B_n fingerprint classification.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from p15_full_bn_prespecht_fixed_plus_scout import (
    block_matrix_mod,
    build_prespecht_terms,
    load_defect_rows,
    mod_rank,
    needed_specht_matrices,
    r_operator_mod,
    r_sigmas,
)
from p15_full_bn_t1_k2_extended_quotient_layer_scout import (
    oriented_pair_lower_boundary,
    oriented_pair_operator,
    part_str,
    specht_rank_for_shape,
    subset_lower,
    subset_operator,
)

PRIME = 1000003
ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ModuleProfile:
    label: str
    target_shape: str
    model: str
    natural_dim: int
    lower_description: str
    lower_column_count: int
    lower_rank: int
    natural_operator_rank: int
    concatenated_lower_image_rank: int
    image_mod_lower_rank: int
    image_intersection_lower_rank: int
    source_lower_image_rank: int
    source_lower_image_outside_lower_rank: int
    specht_dim: int
    specht_block_rank: int
    pass_conditions: Dict[str, bool]


@dataclass(frozen=True)
class RowCheck:
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
class T4Summary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    prime: int
    subset_count: int
    raw_active_terms: int
    raw_sum_abs_coefficients: int
    group_keys: int
    zero_keys_removed: int
    profiles: List[ModuleProfile]
    row_checks: List[RowCheck]
    runtime_seconds: float
    boundary: str
    next_task: str


def profile_from_matrices(
    label: str,
    target_shape: Tuple[int, ...],
    model: str,
    lower_description: str,
    M: np.ndarray,
    B: np.ndarray | None,
    subset_count: int,
    grouped,
    p: int,
    expected_natural_rank: int,
    expected_lower_rank: int,
    expected_image_mod_lower: int,
    expected_source_lower_image: int,
) -> ModuleProfile:
    rank_M = mod_rank(M, p)
    dim, specht_rank = specht_rank_for_shape(target_shape, subset_count, grouped, p)
    if B is None:
        rank_B = 0
        lower_cols = 0
        concat_rank = rank_M
        image_mod_lower = rank_M
        intersection = 0
        source_lower_image = 0
        source_lower_outside = 0
    else:
        rank_B = mod_rank(B, p)
        lower_cols = B.shape[1]
        concat_rank = mod_rank(np.hstack([B, M]), p)
        image_mod_lower = concat_rank - rank_B
        intersection = rank_M + rank_B - concat_rank
        MB = (M @ B) % p
        source_lower_image = mod_rank(MB, p)
        source_lower_outside = mod_rank(np.hstack([B, MB]), p) - rank_B
    checks = {
        "natural_rank_matches_expected": rank_M == expected_natural_rank,
        "lower_rank_matches_expected": rank_B == expected_lower_rank,
        "image_mod_lower_matches_expected": image_mod_lower == expected_image_mod_lower,
        "source_lower_image_rank_matches_expected": source_lower_image == expected_source_lower_image,
        "source_lower_maps_inside_lower": source_lower_outside == 0,
        "specht_rank_matches_quotient_rank": specht_rank == image_mod_lower,
    }
    return ModuleProfile(
        label=label,
        target_shape=part_str(target_shape),
        model=model,
        natural_dim=M.shape[0],
        lower_description=lower_description,
        lower_column_count=lower_cols,
        lower_rank=rank_B,
        natural_operator_rank=rank_M,
        concatenated_lower_image_rank=concat_rank,
        image_mod_lower_rank=image_mod_lower,
        image_intersection_lower_rank=intersection,
        source_lower_image_rank=source_lower_image,
        source_lower_image_outside_lower_rank=source_lower_outside,
        specht_dim=dim,
        specht_block_rank=specht_rank,
        pass_conditions=checks,
    )


def row_checks(n: int, k: int, plus_size: int, plus_masks: List[int], grouped, p: int) -> List[RowCheck]:
    defect_rows = load_defect_rows()
    needed = {sigma for (_dst, _src, sigma) in grouped.keys()} | r_sigmas(n, plus_size, plus_masks)
    shapes = [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
    out: List[RowCheck] = []
    for shape in shapes:
        lambda_minus = part_str(shape)
        mats, dim = needed_specht_matrices(shape, needed, p)
        M = block_matrix_mod(grouped, len(plus_masks), mats, dim, p)
        R = r_operator_mod(n, plus_size, plus_masks, mats, dim, p)
        total_dim = len(plus_masks) * dim
        rank = mod_rank(M, p)
        bound = mod_rank((np.eye(total_dim, dtype=np.int64) + R) % p, p)
        defect = bound - rank
        bucket = "ZERO" if rank == 0 and bound > 0 else ("EXTRA_DEFECT" if defect > 0 else "EQUAL_R_BOUND")
        expected = defect_rows.get((n, k, str(plus_size), lambda_minus))
        listed = expected is not None
        exp_rank = int(expected["rank_mod_p"]) if expected else None
        exp_defect = int(expected["defect_to_bound"]) if expected else None
        exp_bucket = expected["bucket"] if expected else None
        if expected:
            csv_match = (rank == exp_rank and defect == exp_defect and bucket == exp_bucket)
        else:
            csv_match = None if bucket == "EQUAL_R_BOUND" else False
        out.append(RowCheck(
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
    return out


def certify(out_json: Path | None, out_md: Path | None) -> T4Summary:
    t0 = time.time()
    n = 8
    k = 3
    plus_size = 4
    plus_masks, grouped, _perms, raw_active, raw_abs, zero_removed = build_prespecht_terms(n, k, plus_size)
    subset_count = len(plus_masks)
    profiles = [
        profile_from_matrices(
            label="trivial_complement",
            target_shape=(4,),
            model="0-subsets/constants in each 4-point complement fiber",
            lower_description="none",
            M=subset_operator(n, 0, plus_masks, grouped, PRIME),
            B=None,
            subset_count=subset_count,
            grouped=grouped,
            p=PRIME,
            expected_natural_rank=26,
            expected_lower_rank=0,
            expected_image_mod_lower=26,
            expected_source_lower_image=0,
        ),
        profile_from_matrices(
            label="subset_r1_standard",
            target_shape=(3, 1),
            model="1-subsets in each 4-point complement fiber",
            lower_description="fiber constants U_0",
            M=subset_operator(n, 1, plus_masks, grouped, PRIME),
            B=subset_lower(n, 1, 0, plus_masks),
            subset_count=subset_count,
            grouped=grouped,
            p=PRIME,
            expected_natural_rank=104,
            expected_lower_rank=70,
            expected_image_mod_lower=78,
            expected_source_lower_image=26,
        ),
        profile_from_matrices(
            label="subset_r2_two_row",
            target_shape=(2, 2),
            model="2-subsets in each 4-point complement fiber",
            lower_description="fiber point-star subspace U_1",
            M=subset_operator(n, 2, plus_masks, grouped, PRIME),
            B=subset_lower(n, 2, 1, plus_masks),
            subset_count=subset_count,
            grouped=grouped,
            p=PRIME,
            expected_natural_rank=168,
            expected_lower_rank=280,
            expected_image_mod_lower=64,
            expected_source_lower_image=104,
        ),
        profile_from_matrices(
            label="oriented_pair_hook",
            target_shape=(2, 1, 1),
            model="oriented 2-subsets in each 4-point complement fiber",
            lower_description="boundary image of the point module inside oriented pairs",
            M=oriented_pair_operator(n, plus_masks, grouped, PRIME),
            B=oriented_pair_lower_boundary(n, plus_masks, PRIME),
            subset_count=subset_count,
            grouped=grouped,
            p=PRIME,
            expected_natural_rank=168,
            expected_lower_rank=210,
            expected_image_mod_lower=90,
            expected_source_lower_image=78,
        ),
    ]
    rows = row_checks(n, k, plus_size, plus_masks, grouped, PRIME)
    row_pass = all(row.csv_match is not False for row in rows)
    profile_pass = all(all(profile.pass_conditions.values()) for profile in profiles)
    sign_row_pass = any(row.lambda_minus == "1.1.1.1" and row.bucket == "EQUAL_R_BOUND" and row.rank_mod_p == row.r_plus_bound_mod_p == 38 for row in rows)
    status = "PASS" if profile_pass and row_pass and sign_row_pass else "FAIL"
    summary = T4Summary(
        status=status,
        claim=(
            "The T_4^(8,3) one-row fixed-plus family is finitely accounted: all four "
            "extra-defect complement rows are explained by natural quotient layers, and "
            "the sign complement row is equal to its R=+ bound."
        ),
        n=n,
        k=k,
        lambda_plus="4",
        complement_size=4,
        prime=PRIME,
        subset_count=subset_count,
        raw_active_terms=raw_active,
        raw_sum_abs_coefficients=raw_abs,
        group_keys=len(grouped),
        zero_keys_removed=zero_removed,
        profiles=profiles,
        row_checks=rows,
        runtime_seconds=round(time.time() - t0, 3),
        boundary=(
            "Finite n=8,k=3 fixed-plus lambda_plus=4 mechanism certificate only. It accounts for "
            "the T_4^(8,3) one-row family in the n=8 atlas, but it does not prove an all-n theorem "
            "and does not classify the full native B_n fingerprint."
        ),
        next_task=(
            "Keep this package appendix-scoped; remaining fixed-plus families are future finite checks, "
            "and keep this package appendix-scoped; rational T_1/T_2/T_4 audits are separate optional checks."
        ),
    )
    if out_json:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: T4Summary, path: Path) -> None:
    lines = [
        "# P15 T4 k=3 Quotient Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Compressed Operator",
        "",
        "| field | value |",
        "|---|---:|",
        f"| subset count | {summary.subset_count} |",
        f"| raw active terms | {summary.raw_active_terms} |",
        f"| raw sum abs coefficients | {summary.raw_sum_abs_coefficients} |",
        f"| group keys | {summary.group_keys} |",
        f"| zero keys removed | {summary.zero_keys_removed} |",
        "",
        "## Quotient Profiles",
        "",
        "| target | model | natural rank | lower rank | quotient rank | M(lower) rank | outside lower | Specht rank |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for p in summary.profiles:
        lines.append(
            f"| `{p.target_shape}` | {p.model} | {p.natural_operator_rank} | {p.lower_rank} | "
            f"{p.image_mod_lower_rank} | {p.source_lower_image_rank} | "
            f"{p.source_lower_image_outside_lower_rank} | {p.specht_block_rank} |"
        )
    lines += [
        "",
        "## Row Checks",
        "",
        "| lambda_minus | dim S^nu | total dim | R+ bound | rank mod p | defect | bucket | CSV match |",
        "|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in summary.row_checks:
        match = "NA" if row.csv_match is None else str(row.csv_match)
        lines.append(
            f"| `{row.lambda_minus}` | {row.specht_dim} | {row.total_dim} | {row.r_plus_bound_mod_p} | "
            f"{row.rank_mod_p} | {row.defect_to_bound} | {row.bucket} | {match} |"
        )
    lines += [
        "",
        "The sign complement `(1^4)` is not an extra-defect row: its rank and `R=+` bound are both `38`.",
        "",
        "## Checks",
        "",
    ]
    for p in summary.profiles:
        lines += [
            f"### `{p.target_shape}`",
            "",
            f"Lower space: {p.lower_description}.",
            "",
            "```json",
            json.dumps(p.pass_conditions, indent=2),
            "```",
            "",
        ]
    lines += [
        "## Sources",
        "",
        "```text",
        "scripts/p15_full_bn_t4_k3_quotient_certificate.py",
        "certified/P15_FULL_BN_T4_K3_QUOTIENT_CERTIFICATE.json",
        "certified/P15_FULL_BN_T4_K3_QUOTIENT_CERTIFICATE.md",
        "```",
        "",
        "## Boundary",
        "",
        summary.boundary,
        "",
        "Prossima task: " + summary.next_task,
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-json", default=None)
    ap.add_argument("--write-md", default=None)
    args = ap.parse_args()
    summary = certify(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())