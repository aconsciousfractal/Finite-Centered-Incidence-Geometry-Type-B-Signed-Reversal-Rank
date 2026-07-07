#!/usr/bin/env python3
"""
Quotient-layer and sign-twist certificate for T_2^(8,3).

This handles the fixed-plus family lambda_plus=2 in the n=8,k=3 full native
B_n fingerprint atlas.  The compressed operator is

    T_2^(8,3) in Mat_28(Q[S_6]).

The dominant complement shapes are accounted for by natural complement-fiber
modules, and conjugate partner ranks are transferred by a row/column sign-twist
factorization of the complement sign character.

This is a finite n=8,k=3 certificate. It is not an all-n theorem and not a
full native B_n fingerprint classification.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict, deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from p15_full_bn_prespecht_fixed_plus_scout import (
    build_prespecht_terms,
    mod_rank,
    perm_sign,
)
from p15_full_bn_t1_k2_extended_quotient_layer_scout import (
    flag_lower,
    flag_operator,
    oriented_pair_lower_boundary,
    oriented_pair_operator,
    part_str,
    specht_rank_for_shape,
    subset_lower,
    subset_operator,
)

PRIME = 1000003
ROOT = Path(__file__).resolve().parents[1]
PRESPECHT_JSON = ROOT / "certified" / "P15_FULL_BN_PRESPECHT_FIXED_PLUS_SCOUT.json"


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
class RankPair:
    shape: str
    conjugate: str
    rank: int | None
    conjugate_rank: int | None
    r_plus_bound: int | None
    conjugate_r_plus_bound: int | None
    bucket: str | None
    conjugate_bucket: str | None
    rank_match: bool | None
    reading: str


@dataclass(frozen=True)
class SignTwistProfile:
    group_keys: int
    active_cells: int
    assigned_row_nodes: int
    assigned_column_nodes: int
    unassigned_row_nodes: int
    unassigned_column_nodes: int
    row_sign_plus: int
    row_sign_minus: int
    column_sign_plus: int
    column_sign_minus: int
    factorization_bad_terms: int
    consistency_ok: bool
    conjugate_rank_pairs: List[RankPair]
    pass_conditions: Dict[str, bool]


@dataclass(frozen=True)
class T2Summary:
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
    sign_twist: SignTwistProfile
    runtime_seconds: float
    boundary: str
    next_task: str


def parse_part(s: str) -> Tuple[int, ...]:
    if s == "empty":
        return ()
    return tuple(int(x) for x in s.split("."))


def conjugate(part: Tuple[int, ...]) -> Tuple[int, ...]:
    if not part:
        return ()
    return tuple(sum(1 for x in part if x >= j) for j in range(1, max(part) + 1))


def load_t2_rows() -> Dict[str, Dict[str, int | str]]:
    if not PRESPECHT_JSON.exists():
        return {}
    data = json.loads(PRESPECHT_JSON.read_text(encoding="utf-8"))
    for case in data.get("cases", []):
        if case.get("label") == "T_2^(8,3)":
            out: Dict[str, Dict[str, int | str]] = {}
            for row in case.get("rows", []):
                out[row["lambda_minus"]] = {
                    "rank": int(row["rank_mod_p"]),
                    "bound": int(row["r_plus_bound_mod_p"]),
                    "bucket": row["bucket"],
                }
            return out
    return {}


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


def sign_twist_profile(grouped, subset_count: int) -> SignTwistProfile:
    adjacency: Dict[Tuple[str, int], List[Tuple[Tuple[str, int], int]]] = defaultdict(list)
    cells = set()
    for (dst, src, sigma), _coeff in grouped.items():
        sgn = perm_sign(sigma)
        cells.add((dst, src))
        adjacency[("row", dst)].append((("column", src), sgn))
        adjacency[("column", src)].append((("row", dst), sgn))

    assignment: Dict[Tuple[str, int], int] = {}
    consistent = True
    for node in list(adjacency):
        if node in assignment:
            continue
        assignment[node] = 1
        queue = deque([node])
        while queue:
            u = queue.popleft()
            for v, label in adjacency[u]:
                value = assignment[u] * label
                if v in assignment:
                    if assignment[v] != value:
                        consistent = False
                else:
                    assignment[v] = value
                    queue.append(v)

    row_sign = {i: assignment.get(("row", i), 1) for i in range(subset_count)}
    column_sign = {i: assignment.get(("column", i), 1) for i in range(subset_count)}
    bad_terms = 0
    for (dst, src, sigma), _coeff in grouped.items():
        if perm_sign(sigma) != row_sign[dst] * column_sign[src]:
            bad_terms += 1

    ranks = load_t2_rows()
    selected = ["6", "5.1", "4.2", "4.1.1", "3.3", "3.2.1"]
    pairs: List[RankPair] = []
    for shape in selected:
        conj = part_str(conjugate(parse_part(shape)))
        row = ranks.get(shape, {})
        crow = ranks.get(conj, {})
        rank = row.get("rank") if row else None
        crank = crow.get("rank") if crow else None
        bucket = row.get("bucket") if row else None
        cbucket = crow.get("bucket") if crow else None
        if shape == conj:
            reading = "self-conjugate dominant shape; quotient profile explains the rank directly"
        elif shape == "6":
            reading = "rank transfers to the sign row, whose R=+ bound is also 12, so it is equal rather than extra-defective"
        else:
            reading = "partner rank follows from the dominant-side quotient profile by sign-twist"
        pairs.append(RankPair(
            shape=shape,
            conjugate=conj,
            rank=int(rank) if rank is not None else None,
            conjugate_rank=int(crank) if crank is not None else None,
            r_plus_bound=int(row.get("bound")) if row else None,
            conjugate_r_plus_bound=int(crow.get("bound")) if crow else None,
            bucket=str(bucket) if bucket is not None else None,
            conjugate_bucket=str(cbucket) if cbucket is not None else None,
            rank_match=(rank == crank) if rank is not None and crank is not None else None,
            reading=reading,
        ))

    assigned_rows = sum(1 for i in range(subset_count) if ("row", i) in assignment)
    assigned_cols = sum(1 for i in range(subset_count) if ("column", i) in assignment)
    checks = {
        "factorization_consistent": consistent,
        "factorization_bad_terms_zero": bad_terms == 0,
        "all_conjugate_rank_pairs_match": all(p.rank_match is not False for p in pairs),
    }
    return SignTwistProfile(
        group_keys=len(grouped),
        active_cells=len(cells),
        assigned_row_nodes=assigned_rows,
        assigned_column_nodes=assigned_cols,
        unassigned_row_nodes=subset_count - assigned_rows,
        unassigned_column_nodes=subset_count - assigned_cols,
        row_sign_plus=sum(1 for v in row_sign.values() if v == 1),
        row_sign_minus=sum(1 for v in row_sign.values() if v == -1),
        column_sign_plus=sum(1 for v in column_sign.values() if v == 1),
        column_sign_minus=sum(1 for v in column_sign.values() if v == -1),
        factorization_bad_terms=bad_terms,
        consistency_ok=consistent,
        conjugate_rank_pairs=pairs,
        pass_conditions=checks,
    )


def certify(out_json: Path | None, out_md: Path | None) -> T2Summary:
    t0 = time.time()
    n = 8
    k = 3
    plus_size = 2
    plus_masks, grouped, _perms, _raw_active, _raw_abs, _zero_removed = build_prespecht_terms(n, k, plus_size)
    subset_count = len(plus_masks)
    profiles: List[ModuleProfile] = []

    profiles.append(profile_from_matrices(
        label="trivial_complement",
        target_shape=(6,),
        model="0-subsets/constants in each 6-point complement fiber",
        lower_description="none",
        M=subset_operator(n, 0, plus_masks, grouped, PRIME),
        B=None,
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=12,
        expected_lower_rank=0,
        expected_image_mod_lower=12,
        expected_source_lower_image=0,
    ))
    profiles.append(profile_from_matrices(
        label="subset_r1_standard",
        target_shape=(5, 1),
        model="1-subsets in each 6-point complement fiber",
        lower_description="fiber constants U_0",
        M=subset_operator(n, 1, plus_masks, grouped, PRIME),
        B=subset_lower(n, 1, 0, plus_masks),
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=48,
        expected_lower_rank=28,
        expected_image_mod_lower=36,
        expected_source_lower_image=12,
    ))
    profiles.append(profile_from_matrices(
        label="subset_r2_two_row",
        target_shape=(4, 2),
        model="2-subsets in each 6-point complement fiber",
        lower_description="fiber point-star subspace U_1",
        M=subset_operator(n, 2, plus_masks, grouped, PRIME),
        B=subset_lower(n, 2, 1, plus_masks),
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=108,
        expected_lower_rank=168,
        expected_image_mod_lower=60,
        expected_source_lower_image=48,
    ))
    profiles.append(profile_from_matrices(
        label="oriented_pair_hook",
        target_shape=(4, 1, 1),
        model="oriented 2-subsets in each 6-point complement fiber",
        lower_description="boundary image of the point module inside oriented pairs",
        M=oriented_pair_operator(n, plus_masks, grouped, PRIME),
        B=oriented_pair_lower_boundary(n, plus_masks, PRIME),
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=84,
        expected_lower_rank=140,
        expected_image_mod_lower=48,
        expected_source_lower_image=36,
    ))
    profiles.append(profile_from_matrices(
        label="subset_r3_three_row",
        target_shape=(3, 3),
        model="3-subsets in each 6-point complement fiber",
        lower_description="fiber pair-incidence subspace U_2",
        M=subset_operator(n, 3, plus_masks, grouped, PRIME),
        B=subset_lower(n, 3, 2, plus_masks),
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=144,
        expected_lower_rank=420,
        expected_image_mod_lower=36,
        expected_source_lower_image=108,
    ))
    profiles.append(profile_from_matrices(
        label="pair_point_flag_mixed",
        target_shape=(3, 2, 1),
        model="nonoriented pair-point flags (pair P, point z notin P) in each complement fiber",
        lower_description="span of point, pair, triple, and oriented-pair hook incidence submodules",
        M=flag_operator(n, plus_masks, grouped, PRIME),
        B=flag_lower(n, plus_masks, PRIME),
        subset_count=subset_count,
        grouped=grouped,
        p=PRIME,
        expected_natural_rank=384,
        expected_lower_rank=1232,
        expected_image_mod_lower=96,
        expected_source_lower_image=288,
    ))

    twist = sign_twist_profile(grouped, subset_count)
    status = "PASS" if all(all(p.pass_conditions.values()) for p in profiles) and all(twist.pass_conditions.values()) else "FAIL"
    summary = T2Summary(
        status=status,
        claim=(
            "The T_2^(8,3) fixed-plus family is finitely accounted: six dominant-side "
            "complement rows are explained by natural quotient layers, and the nontrivial "
            "conjugate partner ranks transfer by complement sign-twist row/column factorization."
        ),
        n=n,
        k=k,
        lambda_plus="2",
        complement_size=6,
        prime=PRIME,
        subset_count=subset_count,
        group_keys=len(grouped),
        profiles=profiles,
        sign_twist=twist,
        runtime_seconds=round(time.time() - t0, 3),
        boundary=(
            "Finite n=8,k=3 fixed-plus lambda_plus=2 mechanism certificate only. It accounts for "
            "the T_2^(8,3) compact family in the n=8 atlas, but it does not prove an all-n theorem "
            "and does not classify the full native B_n fingerprint."
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


def write_markdown(summary: T2Summary, path: Path) -> None:
    lines = [
        "# P15 T2 k=3 Quotient and Sign-Twist Certificate",
        "",
        f"Status: **{summary.status}**",
        "",
        "## Claim",
        "",
        "```text",
        summary.claim,
        "```",
        "",
        "## Quotient Profiles",
        "",
        "| target | model | natural rank | lower rank | quotient rank | M(lower) rank | outside lower | Specht rank | bucket consequence |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    bucket_notes = {
        "6": "extra-defect row; sign partner is equal-bound",
        "5.1": "extra-defect dominant row",
        "4.2": "extra-defect dominant row",
        "4.1.1": "extra-defect dominant row",
        "3.3": "extra-defect dominant row",
        "3.2.1": "self-conjugate extra-defect row",
    }
    for p in summary.profiles:
        lines.append(
            f"| `{p.target_shape}` | {p.model} | {p.natural_operator_rank} | {p.lower_rank} | "
            f"{p.image_mod_lower_rank} | {p.source_lower_image_rank} | "
            f"{p.source_lower_image_outside_lower_rank} | {p.specht_block_rank} | {bucket_notes[p.target_shape]} |"
        )
    lines += [
        "",
        "## Sign-Twist Transfer",
        "",
        "The complement sign character factors on every nonzero compressed term as a row sign times a column sign on the 28 plus fibers. Four row fibers and four column fibers are isolated in this compressed operator; assigning them sign `+1` is harmless because no nonzero term uses them.",
        "",
        "| field | value |",
        "|---|---:|",
        f"| group keys | {summary.sign_twist.group_keys} |",
        f"| active cells | {summary.sign_twist.active_cells} |",
        f"| assigned row nodes | {summary.sign_twist.assigned_row_nodes} |",
        f"| assigned column nodes | {summary.sign_twist.assigned_column_nodes} |",
        f"| unassigned row nodes | {summary.sign_twist.unassigned_row_nodes} |",
        f"| unassigned column nodes | {summary.sign_twist.unassigned_column_nodes} |",
        f"| bad factorization terms | {summary.sign_twist.factorization_bad_terms} |",
        "",
        "| shape | conjugate | rank | bound | bucket | conjugate rank | conjugate bound | conjugate bucket | reading |",
        "|---|---|---:|---:|---|---:|---:|---|---|",
    ]
    for p in summary.sign_twist.conjugate_rank_pairs:
        lines.append(
            f"| `{p.shape}` | `{p.conjugate}` | {p.rank} | {p.r_plus_bound} | {p.bucket} | "
            f"{p.conjugate_rank} | {p.conjugate_r_plus_bound} | {p.conjugate_bucket} | {p.reading} |"
        )
    lines += [
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
        "### sign-twist",
        "",
        "```json",
        json.dumps(summary.sign_twist.pass_conditions, indent=2),
        "```",
        "",
        "## Sources",
        "",
        "```text",
        "scripts/p15_full_bn_t2_k3_quotient_and_sign_twist_certificate.py",
        "certified/P15_FULL_BN_T2_K3_QUOTIENT_AND_SIGN_TWIST_CERTIFICATE.json",
        "certified/P15_FULL_BN_T2_K3_QUOTIENT_AND_SIGN_TWIST_CERTIFICATE.md",
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