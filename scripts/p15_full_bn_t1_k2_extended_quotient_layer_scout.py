#!/usr/bin/env python3
"""
Extended quotient-layer scout for T_1^(8,2) in the P15 n=8 defect atlas.

This is the second fixed-plus quotient scout for the lambda_plus=1 family.
It extends the first two quotient profiles

    (6,1) from the 1-subset module;
    (5,2) from the 2-subset module;

to the next three dominant-side extra-defect shapes

    (5,1,1) from oriented pairs;
    (4,3) from 3-subsets;
    (4,2,1) from nonoriented pair-point flags.

The output is a finite n=8,k=2 mechanism scout. It is not an all-n theorem and
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

from p15_full_bn_prespecht_fixed_plus_scout import (
    block_matrix_mod,
    build_prespecht_terms,
    mod_rank,
    needed_specht_matrices,
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
class ExtendedScoutSummary:
    status: str
    claim: str
    n: int
    k: int
    lambda_plus: str
    complement_size: int
    prime: int
    subset_count: int
    group_keys: int
    profiles: List[ModuleProfile]
    runtime_seconds: float
    boundary: str
    next_task: str


def part_str(part: Tuple[int, ...]) -> str:
    return ".".join(str(x) for x in part) if part else "empty"


def complement_points(n: int, plus_mask: int) -> List[int]:
    return [x for x in range(n) if not ((plus_mask >> x) & 1)]


def profile_from_matrices(
    label: str,
    target_shape: Tuple[int, ...],
    model: str,
    lower_description: str,
    M: np.ndarray,
    B: np.ndarray,
    subset_count: int,
    grouped,
    p: int,
    expected_natural_rank: int,
    expected_lower_rank: int,
    expected_image_mod_lower: int,
    expected_source_lower_image: int,
) -> ModuleProfile:
    rank_M = mod_rank(M, p)
    rank_B = mod_rank(B, p)
    concat_rank = mod_rank(np.hstack([B, M]), p)
    image_mod_lower = concat_rank - rank_B
    intersection = rank_M + rank_B - concat_rank
    MB = (M @ B) % p
    source_lower_image = mod_rank(MB, p)
    source_lower_outside = mod_rank(np.hstack([B, MB]), p) - rank_B
    specht_dim, specht_rank = specht_rank_for_shape(target_shape, subset_count, grouped, p)
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
        lower_column_count=B.shape[1],
        lower_rank=rank_B,
        natural_operator_rank=rank_M,
        concatenated_lower_image_rank=concat_rank,
        image_mod_lower_rank=image_mod_lower,
        image_intersection_lower_rank=intersection,
        source_lower_image_rank=source_lower_image,
        source_lower_image_outside_lower_rank=source_lower_outside,
        specht_dim=specht_dim,
        specht_block_rank=specht_rank,
        pass_conditions=checks,
    )


def specht_rank_for_shape(shape: Tuple[int, ...], subset_count: int, grouped, p: int) -> Tuple[int, int]:
    sigmas = {sigma for (_dst, _src, sigma) in grouped}
    mats, dim = needed_specht_matrices(shape, sigmas, p)
    block = block_matrix_mod(grouped, subset_count, mats, dim, p)
    return dim, mod_rank(block, p)


def subset_basis(n: int, r: int, plus_masks: List[int]):
    basis: List[Tuple[int, Tuple[int, ...]]] = []
    index: Dict[Tuple[int, Tuple[int, ...]], int] = {}
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for sub in itertools.combinations(comp, r):
            key = (plus_i, tuple(sub))
            index[key] = len(basis)
            basis.append(key)
    return basis, index


def subset_operator(n: int, r: int, plus_masks: List[int], grouped, p: int) -> np.ndarray:
    basis, index = subset_basis(n, r, plus_masks)
    M = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for (dst, src, sigma), coeff in grouped.items():
        src_comp = complement_points(n, plus_masks[src])
        dst_comp = complement_points(n, plus_masks[dst])
        for sub in itertools.combinations(src_comp, r):
            source_ranks = [src_comp.index(x) for x in sub]
            image = tuple(sorted(dst_comp[sigma[t]] for t in source_ranks))
            M[index[(dst, image)], index[(src, tuple(sub))]] = (
                M[index[(dst, image)], index[(src, tuple(sub))]] + coeff
            ) % p
    return M


def subset_lower(n: int, r: int, lower_size: int, plus_masks: List[int]) -> np.ndarray:
    basis, index = subset_basis(n, r, plus_masks)
    columns: List[np.ndarray] = []
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for low in itertools.combinations(comp, lower_size):
            v = np.zeros((len(basis),), dtype=np.int64)
            low_set = set(low)
            for sub in itertools.combinations(comp, r):
                if low_set.issubset(sub):
                    v[index[(plus_i, tuple(sub))]] = 1
            columns.append(v)
    return np.stack(columns, axis=1)


def oriented_pair_basis(n: int, plus_masks: List[int]):
    basis: List[Tuple[int, Tuple[int, int]]] = []
    index: Dict[Tuple[int, Tuple[int, int]], int] = {}
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for pair in itertools.combinations(comp, 2):
            key = (plus_i, tuple(pair))
            index[key] = len(basis)
            basis.append(key)
    return basis, index


def oriented_pair_operator(n: int, plus_masks: List[int], grouped, p: int) -> np.ndarray:
    basis, index = oriented_pair_basis(n, plus_masks)
    M = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for (dst, src, sigma), coeff in grouped.items():
        src_comp = complement_points(n, plus_masks[src])
        dst_comp = complement_points(n, plus_masks[dst])
        for pair in itertools.combinations(src_comp, 2):
            ranks = [src_comp.index(x) for x in pair]
            values = [dst_comp[sigma[t]] for t in ranks]
            if values[0] < values[1]:
                image = (values[0], values[1])
                sign = 1
            else:
                image = (values[1], values[0])
                sign = -1
            M[index[(dst, image)], index[(src, tuple(pair))]] = (
                M[index[(dst, image)], index[(src, tuple(pair))]] + coeff * sign
            ) % p
    return M


def oriented_pair_lower_boundary(n: int, plus_masks: List[int], p: int) -> np.ndarray:
    basis, index = oriented_pair_basis(n, plus_masks)
    columns: List[np.ndarray] = []
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for x in comp:
            v = np.zeros((len(basis),), dtype=np.int64)
            for y in comp:
                if y == x:
                    continue
                if x < y:
                    pair = (x, y)
                    sign = 1
                else:
                    pair = (y, x)
                    sign = -1
                v[index[(plus_i, pair)]] = (v[index[(plus_i, pair)]] + sign) % p
            columns.append(v)
    return np.stack(columns, axis=1)


def flag_basis(n: int, plus_masks: List[int]):
    basis: List[Tuple[int, Tuple[int, int], int]] = []
    index: Dict[Tuple[int, Tuple[int, int], int], int] = {}
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        for pair in itertools.combinations(comp, 2):
            for z in comp:
                if z in pair:
                    continue
                key = (plus_i, tuple(pair), z)
                index[key] = len(basis)
                basis.append(key)
    return basis, index


def flag_operator(n: int, plus_masks: List[int], grouped, p: int) -> np.ndarray:
    basis, index = flag_basis(n, plus_masks)
    M = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for (dst, src, sigma), coeff in grouped.items():
        src_comp = complement_points(n, plus_masks[src])
        dst_comp = complement_points(n, plus_masks[dst])
        for pair in itertools.combinations(src_comp, 2):
            pair_ranks = [src_comp.index(x) for x in pair]
            image_pair = tuple(sorted(dst_comp[sigma[t]] for t in pair_ranks))
            for z in src_comp:
                if z in pair:
                    continue
                image_z = dst_comp[sigma[src_comp.index(z)]]
                M[index[(dst, image_pair, image_z)], index[(src, tuple(pair), z)]] = (
                    M[index[(dst, image_pair, image_z)], index[(src, tuple(pair), z)]] + coeff
                ) % p
    return M


def flag_lower(n: int, plus_masks: List[int], p: int) -> np.ndarray:
    basis, index = flag_basis(n, plus_masks)
    columns: List[np.ndarray] = []
    for plus_i, plus_mask in enumerate(plus_masks):
        comp = complement_points(n, plus_mask)
        # Point incidence: fix z and sum over all pairs disjoint from z.
        for z in comp:
            v = np.zeros((len(basis),), dtype=np.int64)
            for pair in itertools.combinations([x for x in comp if x != z], 2):
                v[index[(plus_i, tuple(pair), z)]] = 1
            columns.append(v)
        # Pair incidence: fix a pair and sum over all outside points.
        for pair in itertools.combinations(comp, 2):
            v = np.zeros((len(basis),), dtype=np.int64)
            for z in comp:
                if z not in pair:
                    v[index[(plus_i, tuple(pair), z)]] = 1
            columns.append(v)
        # Triple incidence: fix a triple and let each point play the singleton role.
        for triple in itertools.combinations(comp, 3):
            v = np.zeros((len(basis),), dtype=np.int64)
            for z in triple:
                pair = tuple(x for x in triple if x != z)
                v[index[(plus_i, pair, z)]] = 1
            columns.append(v)
        # Oriented-pair hook incidence: a wedge b maps to sum_z ({a,z},b)-({b,z},a).
        for a, b in itertools.combinations(comp, 2):
            v = np.zeros((len(basis),), dtype=np.int64)
            for z in comp:
                if z == a or z == b:
                    continue
                pair_a = tuple(sorted((a, z)))
                pair_b = tuple(sorted((b, z)))
                v[index[(plus_i, pair_a, b)]] = (v[index[(plus_i, pair_a, b)]] + 1) % p
                v[index[(plus_i, pair_b, a)]] = (v[index[(plus_i, pair_b, a)]] - 1) % p
            columns.append(v)
    return np.stack(columns, axis=1)


def certify(out_json: Path | None, out_md: Path | None) -> ExtendedScoutSummary:
    t0 = time.time()
    n = 8
    k = 2
    plus_masks, grouped, _perms, _raw_active, _raw_abs, _zero_removed = build_prespecht_terms(n, k, 1)
    profiles: List[ModuleProfile] = []

    profiles.append(profile_from_matrices(
        label="subset_r1_standard",
        target_shape=(6, 1),
        model="1-subsets in each 7-point complement fiber",
        lower_description="fiber constants U_0",
        M=subset_operator(n, 1, plus_masks, grouped, PRIME),
        B=subset_lower(n, 1, 0, plus_masks),
        subset_count=len(plus_masks),
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=12,
        expected_lower_rank=8,
        expected_image_mod_lower=12,
        expected_source_lower_image=0,
    ))
    profiles.append(profile_from_matrices(
        label="subset_r2_two_row",
        target_shape=(5, 2),
        model="2-subsets in each 7-point complement fiber",
        lower_description="fiber point-star subspace U_1",
        M=subset_operator(n, 2, plus_masks, grouped, PRIME),
        B=subset_lower(n, 2, 1, plus_masks),
        subset_count=len(plus_masks),
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=48,
        expected_lower_rank=56,
        expected_image_mod_lower=36,
        expected_source_lower_image=12,
    ))
    profiles.append(profile_from_matrices(
        label="oriented_pair_hook",
        target_shape=(5, 1, 1),
        model="oriented 2-subsets in each 7-point complement fiber",
        lower_description="boundary image of the point module inside oriented pairs",
        M=oriented_pair_operator(n, plus_masks, grouped, PRIME),
        B=oriented_pair_lower_boundary(n, plus_masks, PRIME),
        subset_count=len(plus_masks),
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=60,
        expected_lower_rank=48,
        expected_image_mod_lower=48,
        expected_source_lower_image=12,
    ))
    profiles.append(profile_from_matrices(
        label="subset_r3_two_row",
        target_shape=(4, 3),
        model="3-subsets in each 7-point complement fiber",
        lower_description="fiber pair-incidence subspace U_2",
        M=subset_operator(n, 3, plus_masks, grouped, PRIME),
        B=subset_lower(n, 3, 2, plus_masks),
        subset_count=len(plus_masks),
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=84,
        expected_lower_rank=168,
        expected_image_mod_lower=36,
        expected_source_lower_image=48,
    ))
    profiles.append(profile_from_matrices(
        label="pair_point_flag_mixed",
        target_shape=(4, 2, 1),
        model="nonoriented pair-point flags (pair P, point z notin P) in each complement fiber",
        lower_description="span of point, pair, triple, and oriented-pair hook incidence submodules",
        M=flag_operator(n, plus_masks, grouped, PRIME),
        B=flag_lower(n, plus_masks, PRIME),
        subset_count=len(plus_masks),
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=288,
        expected_lower_rank=560,
        expected_image_mod_lower=108,
        expected_source_lower_image=180,
    ))

    status = "PASS" if all(all(p.pass_conditions.values()) for p in profiles) else "FAIL"
    summary = ExtendedScoutSummary(
        status=status,
        claim=(
            "The dominant-side nonzero lambda_plus=1 extra-defect rows of T_1^(8,2) "
            "for shapes (6,1), (5,2), (5,1,1), (4,3), and (4,2,1) are accounted "
            "for by natural quotient layers. In each case the image modulo the declared "
            "lower incidence subspace has the same modular rank as the corresponding "
            "Specht block."
        ),
        n=n,
        k=k,
        lambda_plus="1",
        complement_size=7,
        prime=PRIME,
        subset_count=len(plus_masks),
        group_keys=len(grouped),
        profiles=profiles,
        runtime_seconds=round(time.time() - t0, 3),
        boundary=(
            "Finite n=8,k=2 modular quotient-layer scout only. This accounts for five "
            "dominant-side nonzero lambda_plus=1 extra-defect shapes; later sign-twist "
            "and residual certificates complete the finite T_1^(8,2) accounting. It "
            "does not prove an all-n full-fingerprint theorem."
        ),
        next_task=(
            "Attack the next compact one-row fixed-plus family T_4^(8,3) (lambda_plus=4, n=8,k=3), "
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


def write_markdown(summary: ExtendedScoutSummary, path: Path) -> None:
    lines = [
        "# P15 T1 k=2 Extended Quotient-Layer Scout",
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
        "| target | model | natural rank | lower rank | image mod lower | image cap lower | M(lower) rank | outside lower | Specht rank |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for p in summary.profiles:
        lines.append(
            f"| `{p.target_shape}` | {p.model} | {p.natural_operator_rank} | {p.lower_rank} | "
            f"{p.image_mod_lower_rank} | {p.image_intersection_lower_rank} | "
            f"{p.source_lower_image_rank} | {p.source_lower_image_outside_lower_rank} | {p.specht_block_rank} |"
        )
    lines += [
        "",
        "## Reading",
        "",
        "Each row builds a natural complement-fiber module for `T_1^(8,2)`, declares a lower incidence subspace, checks that the operator maps the lower source space back into the lower target space, and compares the quotient image rank with the direct Specht block rank.",
        "",
        "The five quotient identities banked here are:",
        "",
        "```text",
        "(6,1):     1-subsets modulo fiber constants -> rank 12.",
        "(5,2):     2-subsets modulo point-stars -> rank 36.",
        "(5,1,1):   oriented 2-subsets modulo point-boundaries -> rank 48.",
        "(4,3):     3-subsets modulo pair-incidence -> rank 36.",
        "(4,2,1):   pair-point flags modulo point/pair/triple/oriented-hook incidences -> rank 108.",
        "```",
        "",
        "For `(4,2,1)`, the flag module has natural rank 288. Its declared lower incidence subspace has rank 560, the operator maps this lower space into itself with image rank 180, and the quotient image has rank 108, exactly the Specht block rank for `(4,2,1)`.",
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
        "scripts/p15_full_bn_t1_k2_extended_quotient_layer_scout.py",
        "artifacts/certified/P15_FULL_BN_T1_K2_EXTENDED_QUOTIENT_LAYER_SCOUT.json",
        "artifacts/certified/P15_FULL_BN_T1_K2_EXTENDED_QUOTIENT_LAYER_SCOUT.md",
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
