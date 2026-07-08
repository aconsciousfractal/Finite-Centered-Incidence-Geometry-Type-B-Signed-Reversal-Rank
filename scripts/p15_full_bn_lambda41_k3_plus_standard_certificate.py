#!/usr/bin/env python3
"""
Plus-standard quotient certificate for lambda_plus=4.1 at n=8,k=3.

This handles the full fixed-plus family

    lambda_plus=(4,1), lambda_minus |- 3

in the n=8,k=3 native B_n defect atlas.  The plus Specht module S^(4,1)
is realized as the quotient of the 5-point permutation module by fiber
constants.  For each complement Specht shape nu |- 3, the layer operator
preserves that lower constant subspace, and the quotient image rank equals the
direct Specht block rank.

This is a finite n=8,k=3 mechanism certificate. It is not an all-n theorem and
not a full native B_n fingerprint classification.
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

from p15_full_bn_zero_block_certificate import masks_of_size, restricted_perm_data, rho_perm, sign_poly_target_coeff
from p15_full_bn_prespecht_fixed_plus_scout import load_defect_rows, mod_rank, needed_specht_matrices, part_str

PRIME = 1000003
ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class PlusStandardProfile:
    lambda_minus: str
    model: str
    natural_dim: int
    lower_description: str
    lower_rank: int
    natural_operator_rank: int
    concatenated_lower_image_rank: int
    image_mod_lower_rank: int
    image_intersection_lower_rank: int
    source_lower_image_rank: int
    source_lower_image_outside_lower_rank: int
    specht_plus_dim: int
    specht_minus_dim: int
    direct_specht_dim: int
    direct_specht_rank: int
    r_plus_bound_mod_p: int
    defect_to_bound: int
    bucket: str
    listed_in_defect_csv: bool
    expected_rank_mod_p: int | None
    expected_defect_to_bound: int | None
    expected_bucket: str | None
    csv_match: bool | None
    pass_conditions: Dict[str, bool]


@dataclass(frozen=True)
class Lambda41Summary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    prime: int
    subset_count: int
    group_keys: int
    raw_active_terms: int
    raw_sum_abs_coefficients: int
    zero_keys_removed: int
    profiles: List[PlusStandardProfile]
    runtime_seconds: float
    boundary: str
    next_task: str


def complement_points(n: int, plus_mask: int) -> List[int]:
    return [x for x in range(n) if (plus_mask >> x) & 1]


def build_bipartite_grouped(n: int, k: int, plus_size: int):
    target = (k, k)
    plus_masks = masks_of_size(n, plus_size)
    subset_index = {mask: i for i, mask in enumerate(plus_masks)}
    full_mask = (1 << n) - 1
    grouped: Dict[Tuple[int, int, Tuple[int, ...], Tuple[int, ...]], int] = {}
    raw_active = 0
    raw_abs = 0
    zero_removed = 0
    restrict_cache = {}
    for pi0 in itertools.permutations(range(n)):
        pi = tuple(pi0)
        for src_i, A_mask in enumerate(plus_masks):
            B_mask = full_mask ^ A_mask
            coeff = sign_poly_target_coeff(pi, B_mask, target)
            if coeff == 0:
                continue
            raw_active += 1
            raw_abs += abs(coeff)
            data_a = restrict_cache.get((pi, A_mask))
            if data_a is None:
                data_a = restricted_perm_data(pi, A_mask, n)
                restrict_cache[(pi, A_mask)] = data_a
            data_b = restrict_cache.get((pi, B_mask))
            if data_b is None:
                data_b = restricted_perm_data(pi, B_mask, n)
                restrict_cache[(pi, B_mask)] = data_b
            Ap_mask, sigma_a = data_a
            _Bp_mask, sigma_b = data_b
            dst_i = subset_index[Ap_mask]
            key = (dst_i, src_i, sigma_a, sigma_b)
            new_val = grouped.get(key, 0) + coeff
            if new_val:
                grouped[key] = new_val
            elif key in grouped:
                del grouped[key]
                zero_removed += 1
    return plus_masks, grouped, raw_active, raw_abs, zero_removed


def r_minus_sigmas(n: int, plus_masks: List[int]) -> set[Tuple[int, ...]]:
    full_mask = (1 << n) - 1
    rho = rho_perm(n)
    out = set()
    for A_mask in plus_masks:
        B_mask = full_mask ^ A_mask
        _Bp_mask, sigma_b = restricted_perm_data(rho, B_mask, n)
        out.add(sigma_b)
    return out


def plus_standard_operator(n: int, minus_shape: Tuple[int, ...], plus_masks: List[int], grouped, p: int):
    sigmas = {sigma_b for (_dst, _src, _sigma_a, sigma_b) in grouped} | r_minus_sigmas(n, plus_masks)
    minus_mats, minus_dim = needed_specht_matrices(minus_shape, sigmas, p)
    plus_size = len(complement_points(n, plus_masks[0]))
    cell = plus_size * minus_dim
    dim = len(plus_masks) * cell
    M = np.zeros((dim, dim), dtype=np.int64)
    for (dst, src, sigma_a, sigma_b), coeff in grouped.items():
        Z = minus_mats[sigma_b]
        for p_src in range(plus_size):
            p_dst = sigma_a[p_src]
            r0 = dst * cell + p_dst * minus_dim
            c0 = src * cell + p_src * minus_dim
            M[r0:r0 + minus_dim, c0:c0 + minus_dim] = (
                M[r0:r0 + minus_dim, c0:c0 + minus_dim] + coeff * Z
            ) % p
    return M, minus_dim


def lower_constants(plus_size: int, subset_count: int, minus_dim: int) -> np.ndarray:
    cell = plus_size * minus_dim
    dim = subset_count * cell
    columns: List[np.ndarray] = []
    for subset_i in range(subset_count):
        for j in range(minus_dim):
            v = np.zeros((dim,), dtype=np.int64)
            for p_i in range(plus_size):
                v[subset_i * cell + p_i * minus_dim + j] = 1
            columns.append(v)
    return np.stack(columns, axis=1)


def direct_specht_block(n: int, plus_shape: Tuple[int, ...], minus_shape: Tuple[int, ...], plus_masks: List[int], grouped, p: int):
    plus_sigmas = {sigma_a for (_dst, _src, sigma_a, _sigma_b) in grouped}
    minus_sigmas = {sigma_b for (_dst, _src, _sigma_a, sigma_b) in grouped} | r_minus_sigmas(n, plus_masks)
    plus_mats, plus_dim = needed_specht_matrices(plus_shape, plus_sigmas, p)
    minus_mats, minus_dim = needed_specht_matrices(minus_shape, minus_sigmas, p)
    cell = plus_dim * minus_dim
    dim = len(plus_masks) * cell
    M = np.zeros((dim, dim), dtype=np.int64)
    for (dst, src, sigma_a, sigma_b), coeff in grouped.items():
        Z = np.kron(plus_mats[sigma_a], minus_mats[sigma_b]) % p
        r0 = dst * cell
        c0 = src * cell
        M[r0:r0 + cell, c0:c0 + cell] = (M[r0:r0 + cell, c0:c0 + cell] + coeff * Z) % p
    return M, plus_dim, minus_dim


def r_plus_bound(n: int, plus_shape: Tuple[int, ...], minus_shape: Tuple[int, ...], plus_masks: List[int], grouped, p: int) -> int:
    plus_sigmas = {sigma_a for (_dst, _src, sigma_a, _sigma_b) in grouped}
    minus_sigmas = {sigma_b for (_dst, _src, _sigma_a, sigma_b) in grouped} | r_minus_sigmas(n, plus_masks)
    full_mask = (1 << n) - 1
    rho = rho_perm(n)
    subset_index = {mask: i for i, mask in enumerate(plus_masks)}
    for A_mask in plus_masks:
        Ap_mask, sigma_a = restricted_perm_data(rho, A_mask, n)
        plus_sigmas.add(sigma_a)
        B_mask = full_mask ^ A_mask
        _Bp_mask, sigma_b = restricted_perm_data(rho, B_mask, n)
        minus_sigmas.add(sigma_b)
    plus_mats, plus_dim = needed_specht_matrices(plus_shape, plus_sigmas, p)
    minus_mats, minus_dim = needed_specht_matrices(minus_shape, minus_sigmas, p)
    cell = plus_dim * minus_dim
    dim = len(plus_masks) * cell
    R = np.zeros((dim, dim), dtype=np.int64)
    for src, A_mask in enumerate(plus_masks):
        B_mask = full_mask ^ A_mask
        Ap_mask, sigma_a = restricted_perm_data(rho, A_mask, n)
        _Bp_mask, sigma_b = restricted_perm_data(rho, B_mask, n)
        dst = subset_index[Ap_mask]
        Z = np.kron(plus_mats[sigma_a], minus_mats[sigma_b]) % p
        r0 = dst * cell
        c0 = src * cell
        R[r0:r0 + cell, c0:c0 + cell] = (R[r0:r0 + cell, c0:c0 + cell] + Z) % p
    return mod_rank((np.eye(dim, dtype=np.int64) + R) % p, p)


def certify(out_json: Path | None, out_md: Path | None) -> Lambda41Summary:
    t0 = time.time()
    n = 8
    k = 3
    plus_shape = (4, 1)
    plus_size = sum(plus_shape)
    minus_shapes = [(3,), (2, 1), (1, 1, 1)]
    plus_masks, grouped, raw_active, raw_abs, zero_removed = build_bipartite_grouped(n, k, plus_size)
    defect_rows = load_defect_rows()
    profiles: List[PlusStandardProfile] = []
    for minus_shape in minus_shapes:
        lambda_minus = part_str(minus_shape)
        M_nat, minus_dim = plus_standard_operator(n, minus_shape, plus_masks, grouped, PRIME)
        B = lower_constants(plus_size, len(plus_masks), minus_dim)
        rank_M = mod_rank(M_nat, PRIME)
        rank_B = mod_rank(B, PRIME)
        concat_rank = mod_rank(np.hstack([B, M_nat]), PRIME)
        image_mod_lower = concat_rank - rank_B
        intersection = rank_M + rank_B - concat_rank
        MB = (M_nat @ B) % PRIME
        source_lower_image = mod_rank(MB, PRIME)
        source_lower_outside = mod_rank(np.hstack([B, MB]), PRIME) - rank_B
        direct, plus_dim, direct_minus_dim = direct_specht_block(n, plus_shape, minus_shape, plus_masks, grouped, PRIME)
        direct_rank = mod_rank(direct, PRIME)
        bound = r_plus_bound(n, plus_shape, minus_shape, plus_masks, grouped, PRIME)
        defect = bound - direct_rank
        bucket = "ZERO" if direct_rank == 0 and bound > 0 else ("EXTRA_DEFECT" if defect > 0 else "EQUAL_R_BOUND")
        expected = defect_rows.get((n, k, "4.1", lambda_minus))
        listed = expected is not None
        exp_rank = int(expected["rank_mod_p"]) if expected else None
        exp_defect = int(expected["defect_to_bound"]) if expected else None
        exp_bucket = expected["bucket"] if expected else None
        csv_match = None if expected is None else (direct_rank == exp_rank and defect == exp_defect and bucket == exp_bucket)
        checks = {
            "lower_rank_is_expected_fiber_constants": rank_B == len(plus_masks) * minus_dim,
            "source_lower_maps_inside_lower": source_lower_outside == 0,
            "quotient_rank_matches_direct_specht_rank": image_mod_lower == direct_rank,
            "direct_minus_dim_matches_natural_minus_dim": direct_minus_dim == minus_dim,
            "csv_match": csv_match is True,
            "extra_defect": bucket == "EXTRA_DEFECT" and defect > 0,
        }
        profiles.append(PlusStandardProfile(
            lambda_minus=lambda_minus,
            model="plus 5-point permutation module modulo fiber constants, tensored with S^nu on the 3-point complement",
            natural_dim=M_nat.shape[0],
            lower_description="fiber constants in the plus 5-point permutation module",
            lower_rank=rank_B,
            natural_operator_rank=rank_M,
            concatenated_lower_image_rank=concat_rank,
            image_mod_lower_rank=image_mod_lower,
            image_intersection_lower_rank=intersection,
            source_lower_image_rank=source_lower_image,
            source_lower_image_outside_lower_rank=source_lower_outside,
            specht_plus_dim=plus_dim,
            specht_minus_dim=minus_dim,
            direct_specht_dim=direct.shape[0],
            direct_specht_rank=direct_rank,
            r_plus_bound_mod_p=bound,
            defect_to_bound=defect,
            bucket=bucket,
            listed_in_defect_csv=listed,
            expected_rank_mod_p=exp_rank,
            expected_defect_to_bound=exp_defect,
            expected_bucket=exp_bucket,
            csv_match=csv_match,
            pass_conditions=checks,
        ))
    status = "PASS" if all(all(p.pass_conditions.values()) for p in profiles) else "FAIL"
    summary = Lambda41Summary(
        status=status,
        claim="The lambda_plus=(4,1), n=8,k=3 fixed-plus family is finitely accounted: all three complement rows are the plus-standard quotient of the 5-point plus permutation module modulo constants, and their ranks match the defect CSV.",
        n=n,
        k=k,
        lambda_plus="4.1",
        complement_size=n - plus_size,
        prime=PRIME,
        subset_count=len(plus_masks),
        group_keys=len(grouped),
        raw_active_terms=raw_active,
        raw_sum_abs_coefficients=raw_abs,
        zero_keys_removed=zero_removed,
        profiles=profiles,
        runtime_seconds=round(time.time() - t0, 3),
        boundary="Finite n=8,k=3 lambda_plus=(4,1) mechanism certificate only; no all-n fixed-plus theorem and no full native B_n fingerprint theorem.",
        next_task="Optional rational direct-rank audits for T_1/T_2/T_4 are separate replay-tier checks; no full native B_n fingerprint classification is promoted.",
    )
    if out_json is not None:
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")
    if out_md is not None:
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(summary, out_md)
    return summary


def write_markdown(summary: Lambda41Summary, path: Path) -> None:
    lines = [
        "# P15 lambda_plus=4.1 k=3 Plus-Standard Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Profiles",
        "",
        "| lambda_minus | natural rank | lower rank | quotient rank | direct Specht rank | R+ bound | defect | bucket | CSV match |",
        "|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for p in summary.profiles:
        lines.append(
            f"| `{p.lambda_minus}` | {p.natural_operator_rank} | {p.lower_rank} | {p.image_mod_lower_rank} | "
            f"{p.direct_specht_rank} | {p.r_plus_bound_mod_p} | {p.defect_to_bound} | {p.bucket} | {p.csv_match} |"
        )
    lines += [
        "",
        "## Mechanism",
        "",
        "For `lambda_plus=(4,1)`, the plus-side Specht module is the standard quotient of the 5-point permutation module by constants. The certificate builds that natural module over each plus subset fiber, checks that the layer operator maps fiber constants back into fiber constants, and then compares the quotient image rank with the direct `S^(4,1) tensor S^nu` Specht block rank.",
        "",
        "## Checks",
        "",
    ]
    for p in summary.profiles:
        lines += [
            f"### `{p.lambda_minus}`",
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
        "scripts/p15_full_bn_lambda41_k3_plus_standard_certificate.py",
        "artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.json",
        "artifacts/certified/P15_FULL_BN_LAMBDA41_K3_PLUS_STANDARD_CERTIFICATE.md",
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
    parser.add_argument("--write-json", default=None)
    parser.add_argument("--write-md", default=None)
    args = parser.parse_args()
    summary = certify(Path(args.write_json) if args.write_json else None, Path(args.write_md) if args.write_md else None)
    print(json.dumps(asdict(summary), indent=2))
    return 0 if summary.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
